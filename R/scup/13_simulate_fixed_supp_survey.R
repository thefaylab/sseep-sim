### created: 01/18/2024
### updated: 04/15/2025

# 04a - SIMULATE STATUS QUO SURVEY ####


## Objective ####
# For a given species and distribution, simulate the status quo NMFS bottom trawl survey.

# Outputs: one survey and respective tow level data for each replicate of a simulated population and abundance

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
library(purrr)
library(data.table)
library(here)
library(dplyr)
library(tidyverse)
library(cubelyr)
suppressPackageStartupMessages(library(tidyverse))
source(here("R/selectivity_fns.R"))
set.seed(131)


### DATA SET UP ####
# Directories
pop.dat   <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")
dist.dat  <- here("data", "rds", "dists")
survdat   <- here("data", "rds", "survdat")

# Parameters
species   <- "scup"
season    <- "fall"
ages      <- 0:7
years     <- 1:15
nsims     <- 1:100
nsurveys  <- 25
years_post <- 6:15 # Define the years affected by the survey change
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))

# Trawl survey configuration
trawl_dim       <- c(2.7, 0.014)
catch_q         <- force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1)
set_den         <- 0.001
min_sets        <- 3
age_sampling    <- "stratified"
age_space_group <- "set"
resample_cells  <- TRUE

#BLOCK 1
# One population only, inside wind areas only

i <- 1

pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

pop_inside <- pop_i

wind <- pop_inside$grid[["AREA_CODE"]] == 1

for (nm in names(pop_inside$grid)) {
  x <- pop_inside$grid[[nm]]
  x[!wind] <- NA
  pop_inside$grid[[nm]] <- x
}

sum(as.vector(pop_inside$grid$AREA_CODE) == 1, na.rm = TRUE)
table(as.vector(pop_i$grid$AREA_CODE), useNA = "ifany") #original grid

survey_inside_i <- sim_survey(
  pop_inside,
  n_sims = nsurveys,
  trawl_dim = trawl_dim,
  q = catch_q,
  set_den = set_den,
  min_sets = min_sets,
  age_sampling = age_sampling,
  age_space_group = age_space_group,
  resample_cells = TRUE
)

fixed_locs <- survey_inside_i$setdet %>%
  filter(year == 1, sim %in% 1:nsurveys) %>%
  distinct(sim, set, strat, x, y, cell, depth, AREA_CODE,
           tow_area, cell_area, strat_cells, strat_area,
           strat_sets, cell_sets) %>%
  mutate(division = 1)

fixed_inside <- fixed_locs %>%
  tidyr::crossing(year = years_post) %>%
  arrange(sim, year, strat, set) %>%
  mutate(set = row_number()) %>%   # important
  select(sim, year, strat, x, y, cell, depth, AREA_CODE, division,
         tow_area, cell_area, strat_cells, strat_area,
         strat_sets, cell_sets, set)


saveRDS(fixed_inside, here(survdat, "scup_fall_pop001_sims_year01_fixed_inside.rds"))

#BLOCK 2
#RUN survey with fixed locss

# Define chunking
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))

chunk_id <- 5 #change from 1 to 5
this_chunk <- chunks[[chunk_id]]


# c) Load fixed tow locations (created here in this script)
#fixed_inside <- readRDS(here("data", "rds", "survdat", "inside_fixed_locs_survey.rds"))
fixed_inside <- readRDS(here("data", "rds", "survdat", "scup_fall_pop001_sims_year01_fixed_inside.rds"))


# d) Loop over populations in this chunk
set.seed(123)
for (i in this_chunk) {
  message(sprintf("Running fixed supp survey for population %03d of chunk %d...", i, chunk_id))

  # Load population abundance-distribution object
  pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

  # Run survey for fixed stations only (no random sampling, no reallocation)
  survey_fixed <- sim_survey(
    sim              = pop_i,
    n_sims           = nsurveys,           # 25 identical supplemental survey replicates
    trawl_dim        = trawl_dim,
    q                = catch_q,
    resample_cells   = FALSE,              # fixed locations only
    custom_sets      = fixed_inside,       # use your inside-wind grid
    age_sampling     = "random",
    age_space_group  = "set"
  )

  # Save result
  saveRDS(survey_fixed, here(survdat, sprintf("%s_%s_%03d_%d_supp_fixed_survey.rds", species, season, i, nsurveys)))
  message(sprintf("Saved fixed supp survey for population %03d.", i))
}



### BLOCK 3: Arrange supp survey + preclusion
### BLOCK 3: Arrange supp fixed survey + preclusion ####

chunk_size <- 100
chunks <- split(nsims, ceiling(nsims / chunk_size))
chunk_id <- 1
this_chunk <- chunks[[chunk_id]]

for (i in this_chunk) {
  message(sprintf("Combining preclusion + fixed supplemental survey for population %03d...", i))

  precl_survey <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_precl_survey.rds", species, season, i, nsurveys)))
  supp_survey  <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_supp_fixed_survey.rds", species, season, i, nsurveys)))

  precl_setdet <- precl_survey %>%
    mutate(survey_type = "standard_precl")

  supp_setdet <- supp_survey$setdet %>%
    filter(year %in% years_post) %>%
    mutate(survey_type = "supplemental_fixed")

  survey_hybrid <- bind_rows(precl_setdet, supp_setdet) %>%
    arrange(sim, year, strat, set)

  saveRDS(survey_hybrid, here(survdat, sprintf("%s_%s_%03d_%d_hybrid_fixed_survey.rds", species, season, i, nsurveys)))

  message(sprintf("Saved hybrid fixed survey for population %03d: %d rows", i, nrow(survey_hybrid)))
}


