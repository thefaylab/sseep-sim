### created: 01/18/2024
### updated: 07/02/2026

#

## Objective ####


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
# One population only, but now using full strata that CONTAIN wind cells
# (not just the wind cells themselves)

i <- 1

pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

pop_inside <- pop_i

# 1. Identify wind cells
wind <- pop_i$grid[["AREA_CODE"]] == 1



# 2. Restratify: give wind cells their own unique stratum code
# sim_survey() has no direct argument for controlling where it samples. Tt allocates tows (via min_sets/set_den) based entirely on the "strat"
# attribute already stored in the grid object.
# Wind cells currently share a stratum code with the surrounding non-wind area, so sampling inside the wind footprint is left to chance.
# Splitting wind cells into their own stratum here makes min_sets/set_den apply directly to them, guaranteeing tow coverage inside the wind area.


strat_layer <- pop_i$grid[["strat"]]
offset <- 10000

new_strat <- strat_layer
wind_idx <- which(as.vector(wind))
new_strat[wind_idx] <- as.vector(strat_layer)[as.vector(wind_idx)] + offset

pop_inside$grid[["strat"]] <- new_strat

sum(is.na(as.vector(wind)))          # how many NAs in wind
sum(is.na(as.vector(strat_layer)))   # how many NAs in strat_layer

# sanity check: confirms restratification worked as intended
# diagonal entries (row == column) = non-wind cells that keep their orifinal stratum code (unchanged)
# off-diagonal entries where column = row + offset = wind cells that were reassigned to a new stratum code, split out of their original stratum
# counts should match so if stratum 1010 had 192 wind cells, they now appear all under new stratum (1010 + offset), with the  remaining non-wind cells still counted under 1010

table(as.vector(strat_layer), as.vector(new_strat), useNA = "ifany")

survey_inside_i <- sim_survey(
  pop_inside,
  n_sims = nsurveys,
  trawl_dim = trawl_dim,
  q = catch_q,
  set_den = set_den,
  min_sets = min_sets,
  age_sampling = age_sampling,
  age_space_group = age_space_group,
  resample_cells = TRUE)



# 3. Identify original strata that contain any wind cells needed so the supplemental survey pool includes the FULL stratum (wind + non-wind cells sharing that original stratum), not just the wind cells alone

strata_with_wind <- unique(as.vector(strat_layer)[as.vector(wind)])
strata_with_wind <- strata_with_wind[!is.na(strata_with_wind)]

fixed_locs_wind_strata <- survey_inside_i$setdet %>%
  filter(year == 1, sim %in% 1:nsurveys) %>%
  filter(strat %in% strata_with_wind | strat %in% (strata_with_wind + offset)) %>%
  distinct(sim, set, strat, x, y, cell, depth, AREA_CODE,
           tow_area, cell_area, strat_cells, strat_area,
           strat_sets, cell_sets) %>%
  mutate(division = 1)

fixed_inside_wind_strata <- fixed_locs_wind_strata %>%
  tidyr::crossing(year = years_post) %>%
  arrange(sim, year, strat, set) %>%
  mutate(set = row_number()) %>%   # important
  select(sim, year, strat, x, y, cell, depth, AREA_CODE, division,
         tow_area, cell_area, strat_cells, strat_area,
         strat_sets, cell_sets, set)

saveRDS(fixed_inside_wind_strata, here(survdat, "scup_fall_pop001_sims_year01_fixed_inside_wind_strata.rds"))



#BLOCK 2
#RUN survey with fixed locss

# Define chunking
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))

chunk_id <- 1 #change from 1 to 5
this_chunk <- chunks[[chunk_id]]


# c) Load fixed tow locations (created here in this script)
#fixed_inside_wind_strata <- readRDS(here("data", "rds", "survdat", "inside_fixed_locs_wind_strata_survey.rds"))
fixed_inside_wind_strata <- readRDS(here("data", "rds", "survdat", "scup_fall_pop001_sims_year01_fixed_inside_wind_strata.rds"))


# d) Loop over populations in this chunk
set.seed(123)
for (i in this_chunk) {
  message(sprintf("Running fixed supp survey inside wind areas for population %03d of chunk %d...", i, chunk_id))

  # Load population abundance-distribution object
  pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

  # Run survey for fixed stations only (no random sampling, no reallocation)
  survey_fixed_wind_strata <- sim_survey(
    sim              = pop_i,
    n_sims           = nsurveys,           # 25 identical supplemental survey replicates
    trawl_dim        = trawl_dim,
    q                = catch_q,
    resample_cells   = FALSE,              # fixed locations only
    custom_sets      = fixed_inside_wind_strata,       # use inside-wind grid
    age_sampling     = "random",
    age_space_group  = "set"
  )

  # Save result
  saveRDS(survey_fixed_wind_strata, here(survdat, sprintf("%s_%s_%03d_%d_supp_fixed_survey_wind_strata.rds", species, season, i, nsurveys)))
  message(sprintf("Saved fixed supp survey for population %03d.", i))
}




#Check

survey_fixed_wind_strata$setdet |>
  distinct(strat) |>
  filter(!(strat %in% strata_with_wind | strat %in% (strata_with_wind + offset))) #returns 0. any stratum showing up here would mean a tow leaked in from a stratum that does not have wind cells

survey_fixed_wind_strata$setdet %>% count(AREA_CODE) #shows area code 1 and 2 but thats ok because the code intentionally was designed for that.

survey_fixed_wind_strata$setdet |>
  filter(AREA_CODE == 1) |>
  count(strat) #23 strata (new strata) that are wind affected with their respective amount of tows per sim/year/pop
#750 / 25 / 10 = 3 exactly the min_sets - the majority of the strata
#For 11650 (1500 rows): 1500 / 25 / 10 = 6 tows per sim.
#For 11690 (1250 rows): 1250 / 25 / 10 = 5 tows per sim.

fixed_inside_wind_strata %>%
  filter(AREA_CODE == 1) %>%
  distinct(strat, strat_sets, strat_area)


table(survey_fixed_wind_strata$setdet$year) #confirms years 6 to 15 only



### BLOCK 3: Arrange supp fixed survey + preclusion ####
#survdat_supp <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_supp_fixed_survey_wind_strata.rds", species, season, .x))))

chunk_size <- 100
chunks <- split(nsims, ceiling(nsims / chunk_size))
chunk_id <- 1
this_chunk <- chunks[[chunk_id]]

for (i in this_chunk) {
  message(sprintf("Combining preclusion + fixed supplemental survey inside wind strata for population %03d...", i))

  precl_survey <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_precl_survey.rds", species, season, i, nsurveys)))
  supp_wa_survey  <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_supp_fixed_survey_wind_strata.rds", species, season, i, nsurveys)))

  precl_setdet <- precl_survey %>%
    mutate(survey_type = "standard_precl")

  supp_wa_setdet <- supp_wa_survey$setdet %>%
    filter(year %in% years_post) %>%
    mutate(survey_type = "supplemental_fixed_wind_area")

  survey_combined <- bind_rows(precl_setdet, supp_wa_setdet) %>%
    arrange(sim, year, strat, set)

  saveRDS(survey_combined, here(survdat, sprintf("%s_%s_%03d_%d_supp_wa+precl_survey.rds", species, season, i, nsurveys)))

  message(sprintf("Saved supp_wa + precl survey for population %03d: %d rows", i, nrow(survey_combined)))
}


