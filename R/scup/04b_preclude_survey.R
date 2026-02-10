### created: 01/18/2024
### updated: 02/07/2025

# 04b - PRECLUDE SURVEY FROM WIND ENERGY AREAS ####


## Objective ####
# For a given species, distribution, and survey, sets are removed from the simulated tow level data if they occurred in a cell indexed as a wind cell representing an overlap with offshore wind energy areas.

# Outputs: one survey and respective tow level data occurring outside of offshore wind energy areas only for each replicate of a simulated population and abundance
#

### PACKAGES ####
library(tidyverse)
library(here)
library(data.table)


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


# Set which chunk to run (1 to 5)
chunk_id <- 5
this_chunk <- chunks[[chunk_id]]

# Process each chunk
for (i in this_chunk) {
  message(sprintf("Applying preclusion for survey %03d of chunk %d...", i, chunk_id))

  # Load survey result
  surv_i <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_sq_survey.rds", species, season, i, nsurveys)))

  # Apply preclusion (remove AREA_CODE 1 for years >= 6)
  precl_i <- surv_i$setdet |>
    filter(!(year >= 6 & AREA_CODE == 1)) |>
    as.data.table()

  # Save output
  saveRDS(precl_i, here(survdat, sprintf("%s_%s_%03d_%d_precl_survey.rds", species, season, i, nsurveys)))
}



#Compare
survdat_precl_old <- readRDS(here(survdat, "scup_fall_2_sims_25_precl-surv-dat.rds"))
survdat_precl_84 <- readRDS(here(survdat, "scup_fall_084_25_precl_survey.rds"))

survdat_precl_old[[1]]
survdat_precl_84


#Count
s.pre <- survdat_precl_84

# Count sets by sim, year
tow_counts_pre <- s.pre |>
  count(sim, year, name = "n_sets") #can add strat to check

print(tow_counts_pre)


#Average number of tows per year across all simulations
avg_tows_per_year_pre <- s.pre |>
  count(sim, year) |>
  group_by(year) |>
  summarise(
    mean_tows = mean(n),
    sims = n(),
    .groups = "drop"
  )

print(avg_tows_per_year_pre)




#
#
# ### PACKAGES ####
# library(sdmTMB)
# library(SimSurvey)
# suppressPackageStartupMessages(library(tidyverse))
# library(data.table)
# library(here)
# # source(here("R", "sim_pop_fn.R"))
#
#
# ### DATA SET UP ####
# # data locations
# sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
# survdat <- here("data", "rds", "survdat")
#
# ### name of species to be simulated
# species <- "scup"
#
# ### season to be simulated
# season <- "fall"
#
# ### number of simulations
# nsims <- 1:2
#
# ### number of simulations
# nsurveys <- 25
#
# ### ages simulated
# ages <- 0:7
#
# ### years projected
# years <- 1:15
#
#
#
# ### LOAD DATA ####
# # simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
# survdat_sq <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat.rds", sep = "_")))
# #survdat_sq_nw <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat_nw.rds", sep = "_")))
#
# ## PRECLUDE SURVEY ####
# # remove tows/sets that were sampled in cells indexed as a wind cell to represent the preclusion of the survey due to offshore wind energy areas
# # `AREA_CODE == 2` represents outside areas
# survdat_precl <- map(survdat_sq, ~filter(.$setdet, AREA_CODE == 2))
#
# survdat_precl <- map(survdat_sq, ~filter(.$setdet, !(year >= 6 & AREA_CODE == 1)))
# print(table(survdat_precl[[1]]$year, survdat_precl[[1]]$AREA_CODE))  # Check that years 6-15 only have AREA_CODE 2
#
# # Check if AREA_CODE == 1 is removed only for years 6-15
# map(survdat_precl, ~print(table(.x$year, .x$AREA_CODE)))
#
# # Check if AREA_CODE == 1 is removed for years 6-15 across all simulations
# map(survdat_precl, ~table(.x$sim, .x$AREA_CODE, .x$year))
#
#
#
#
# ## SAVE THE DATA ####
# saveRDS(survdat_precl, here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "precl-surv-dat.rds", sep = "_")))
#  #saveRDS(survdat_precl_nw, here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "precl-surv-dat_nw.rds", sep = "_")))
#
