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
map2(trueN_sq, ihat_sq, ~(.y$rel_ihat/.x$rel_N)/.y$rel_ihat)

### TRUE v PRECLUDED ####
map2(trueN_sq, ihat_precl, ~(.y$rel_ihat/.x$rel_N)/.y$rel_ihat)


### TRUE v REALLOCATED ####
map2(trueN_sq, ihat_reall, ~(.y$rel_ihat/.x$rel_N)/.y$rel_ihat)


## ABSOLUTE RELATIVE ERROR ####
### TRUE v STATUS QUO ####
map2(trueN_sq, ihat_sq, ~abs((.y$rel_ihat/.x$rel_N)/.y$rel_ihat))

### TRUE v PRECLUDED ####
map2(trueN_sq, ihat_precl, ~abs((.y$rel_ihat/.x$rel_N)/.y$rel_ihat))

### TRUE v REALLOCATED ####
map2(trueN_sq, ihat_reall, ~abs((.y$rel_ihat/.x$rel_N)/.y$rel_ihat))


