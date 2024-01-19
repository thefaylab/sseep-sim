### created: 01/18/2024
### updated:

# 06a - CALCULATE ACTUAL ERROR ####


## Objective ####
#

# Outputs:

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
# source(here("R", "sim_stratmean_fn.R"))


### DATA SET UP ####
# data locations
sseep.analysis <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis"
survdat <- here("data", "rds", "survdat")
surv.prod <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### LOAD DATA ####
# simulated true abundance created here("R", "05_calculate_rel_abundance.R")
trueN_sq <- readRDS(here(surv.prods, str_c(species, season, "sq_rel-TrueN.rds", sep = "_")))

# simulated status quo abundance index created here("R", "05_calculate_rel_abundance.R")
ihat_sq <- readRDS(here(surv.prods, str_c(species, season, "sq_rel-ihat.rds", sep = "_")))

# simulated precluded abundance index created here("R", "05_calculate_rel_abundance.R")
ihat_precl <- readRDS(here(surv.prods, str_c(species, season, "precl_rel-ihat.rds", sep = "_")))

# simulated reallocated abundance index created here("R", "05_calculate_rel_abundance.R")
ihat_reall <- readRDS(here(surv.prods, str_c(species, season, "reall_rel-ihat.rds", sep = "_")))



## RELATIVE ERROR ####

### TRUE v STATUS QUO ####
true.v.sq_err <- map2(ihat_sq, trueN_sq, ~left_join(.x, .y, by = "year")) |>
  map(~mutate(., rel_err = (.$rel_ihat-.$rel_N)/.$rel_ihat))


### TRUE v PRECLUDED ####
true.v.precl_err <- map2(ihat_precl, trueN_sq, ~left_join(.x, .y, by = "year")) |>
  map(~mutate(., rel_err = (.$rel_ihat-.$rel_N)/.$rel_ihat))


### TRUE v REALLOCATED ####
true.v.reall_err <- map2(ihat_reall, trueN_sq, ~left_join(.x, .y, by = "year")) |>
  map(~mutate(., rel_err = (.$rel_ihat-.$rel_N)/.$rel_ihat))


## ABSOLUTE RELATIVE ERROR ####
### TRUE v STATUS QUO ####
true.v.sq_err <- map(true.v.sq_err, ~mutate(., abs_rel_err = abs(rel_err)))

### TRUE v PRECLUDED ####
true.v.precl_err <- map(true.v.precl_err, ~mutate(., abs_rel_err = abs(rel_err)))

### TRUE v REALLOCATED ####
true.v.reall_err <- map(true.v.reall_err, ~mutate(., abs_rel_err = abs(rel_err)))


## SAVE THE DATA ####
saveRDS(true.v.sq_err, here(perform.metrics, str_c(species, season, "true-v-sq-error.rds", sep = "_")))
saveRDS(true.v.precl_err, here(perform.metrics, str_c(species, season, "true-v-precl-error.rds", sep = "_")))
saveRDS(true.v.reall_err, here(perform.metrics, str_c(species, season, "true-v-reall-error.rds", sep = "_")))


