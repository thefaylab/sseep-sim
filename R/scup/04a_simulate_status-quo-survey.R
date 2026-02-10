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



# Set which chunk to run (1 to 5)
chunk_id <- 1
this_chunk <- chunks[[chunk_id]]

for (i in this_chunk) {
  message(sprintf("Simulating survey for population %03d of chunk %d...", i, chunk_id))

  # Load distributed pop object
  pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

  # Run the survey
  survey_i <- sim_survey(
    pop_i,
    n_sims = nsurveys,
    trawl_dim = trawl_dim,
    q = catch_q,
    set_den = set_den,
    min_sets = min_sets,
    age_sampling = age_sampling,
    age_space_group = age_space_group,
    resample_cells = resample_cells
  )

  # Save survey output
  saveRDS(survey_i, file = here(survdat, sprintf("%s_%s_%03d_25_sq_survey.rds", species, season, i)))
}


#Compare
survdat_sq_old <- readRDS(here(survdat, "scup_fall_2_sims_25_sq-surv-dat.rds"))
survdat_sq_84 <- readRDS(here(survdat, "scup_fall_084_25_sq_survey.rds"))

survdat_sq_old[[1]]$setdet
survdat_sq_84$setdet


#Count
sq <- survdat_sq_84$setdet

# Count sets by sim, year
tow_counts <- sq |>
  count(sim, year, name = "n_sets") #can add strat to check

print(tow_counts)


#Average number of tows per year across all simulations
avg_tows_per_year <- sq |>
  count(sim, year) |>
  group_by(year) |>
  summarise(
    mean_tows = mean(n),
    sims = n(),
    .groups = "drop"
  )

print(avg_tows_per_year)



survdat_sq <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, .x))))


# ### LOAD DATA ####
# # simulated abundance and distributions created here("R", "03_append_distributions.R")
# pop <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "abund-dist.rds", sep = "_")))
# #pop_nw <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "nowind-abund-dist.rds", sep = "_")))
#
#
# ## SIMULATE SURVEY ####
# survdat_sq <- map(pop, ~sim_survey(.$pop,
#                                    n_sims = nsurveys, # one survey per item in population list object
#                                    trawl_dim = trawl_dim,
#                                    q = catch_q,
#                                    set_den = set_den,
#                                    min_sets = min_sets,
#                                    age_sampling = age_sampling,
#                                    age_space_group = age_space_group,
#                                    resample_cells = resample_cells)
# )
#
# # survdat_sq_nw <- map(pop_nw, ~sim_survey(.$pop,
# #                                    n_sims = nsurveys, # one survey per item in population list object
# #                                    trawl_dim = trawl_dim,
# #                                    q = catch_q,
# #                                    set_den = set_den,
# #                                    min_sets = min_sets,
# #                                    age_sampling = age_sampling,
# #                                    age_space_group = age_space_group,
# #                                    resample_cells = resample_cells))
# #
#
# ## SAVE THE DATA ####
# saveRDS(survdat_sq, here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat.rds", sep = "_")))
# #saveRDS(survdat_sq_nw, here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat_nw.rds", sep = "_")))
#
#
#

