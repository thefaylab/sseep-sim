### created: 01/18/2024
### updated:

# 05 - CALCULATE RELATIVE ABUNDANCES ####


## Objective ####
# For a given species, distribution, and survey, calulate the relative true abundnace and the relative abundance index over time.

# Outputs:

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
source(here("R", "sim_stratmean_fn.R"))


### DATA SET UP ####
# data locations
sseep.analysis <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis"
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### LOAD DATA ####
# simulated abundance and distributions created here("R", "03_append_distributions.R")
pop <- readRDS(here(dist.dat, str_c(species, "_abund-dist.rds", sep = "")))

# simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
survdat_sq <- readRDS(here(survdat, str_c(species, season, "sq-survdat.rds", sep = "_")))

# simulated precluded survey data created here("R", "04b_simulate_status-quo-survey.R")
survdat_precl <- readRDS(here(survdat, str_c(species, season, "precl-survdat.rds", sep = "_")))

# simulated reallocated survey data created here("R", "04c_simulate_status-quo-survey.R")
survdat_reall <- readRDS(here(survdat, str_c(species, season, "reall-survdat.rds", sep = "_")))

# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# find the total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))


## True Abundance ####
# calculate the relative true abundance from the simulated population and distribution
trueN_sq <- map(pop, ~as_tibble(.$pop$N) |>
                mutate(age = ages) |>
                pivot_longer(cols = all_of(years),
                             names_to = "year",
                             values_to = "N") |>
                summarise(N = sum(N), .by = "year") |> # calculate the sum of N across ages
                mutate(rel_N = N/mean(N),
                       year = as.integer(year)) # standardize the annual population by the average population size over the projection
              )

## Abundance Index ####
# calculate the abundance index and relative abundance index for each of the scenarios

### Status Quo ####
ihat_sq <- map(survdat_sq, ~as_tibble(.$setdet) |> sim_stratmean( strata_wts = strata_wts, survey_area = survey_area) |>
              mutate(rel_ihat = stratmu/mean(stratmu)))

### Precluded Survey ####
ihat_precl <- map(survdat_precl, ~as_tibble(.) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
                 mutate(rel_ihat = stratmu/mean(stratmu)))

### Reallocated Survey ####
ihat_reall <- map(survdat_reall, ~as_tibble(.x) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
                 mutate(rel_ihat = stratmu/mean(stratmu)))


## SAVE THE DATA ####
saveRDS(trueN_sq, here(surv.prods, str_c(species, season, "sq_rel-TrueN.rds", sep = "_")))
saveRDS(ihat_sq, here(surv.prods, str_c(species, season, "sq_rel-ihat.rds", sep = "_")))
saveRDS(ihat_precl, here(surv.prods, str_c(species, season, "precl_rel-ihat.rds", sep = "_")))
saveRDS(ihat_reall, here(surv.prods, str_c(species, season, "reall_rel-ihat.rds", sep = "_")))
