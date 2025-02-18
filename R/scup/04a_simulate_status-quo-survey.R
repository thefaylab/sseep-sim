### created: 01/18/2024
### updated: 01/21/2025

# 04a - SIMULATE STATUS QUO SURVEY ####


## Objective ####
# For a given species and distribution, simulate the status quo NMFS bottom trawl survey.

# Outputs: one survey and respective tow level data for each replicate of a simulated population and abundance

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
library(purrr)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
# source(here("R", "sim_pop_fn.R"))
set.seed(131)


### DATA SET UP ####
# data locations
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
pop.dat <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")

### name of species to be simulated
species <- "scup"

### season to be simulated
season <- "fall"

### ages simulated
ages <- 0:7

### years projected
years <- 1:15

### number of simulations
nsims <- 1:2

### number of simulations
nsurveys <- 25

### trawl dimensions
trawl_dim <- c(2.7, 0.014)

### catchability from the most recent stock assessment
source(here("R/selectivity_fns.R"))

#catch_q <- sim_logistic(k = 2, x0 = 2.5)
catch_q <- force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1)
### tow density allocation
set_den <- 0.001

### minimum sets to allocate per strata
min_sets <- 3

### age sampling method
age_sampling <- "stratified"

### spatial scale of age sampling
age_space_group <- "set"  # age sampling with with tow and year

### allow resampling of grid cells
resample_cells <- TRUE


### LOAD DATA ####
# simulated abundance and distributions created here("R", "03_append_distributions.R")
pop <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "abund-dist.rds", sep = "_")))
#pop_nw <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "nowind-abund-dist.rds", sep = "_")))


## SIMULATE SURVEY ####
survdat_sq <- map(pop, ~sim_survey(.$pop,
                                   n_sims = nsurveys, # one survey per item in population list object
                                   trawl_dim = trawl_dim,
                                   q = catch_q,
                                   set_den = set_den,
                                   min_sets = min_sets,
                                   age_sampling = age_sampling,
                                   age_space_group = age_space_group,
                                   resample_cells = resample_cells)
)

# survdat_sq_nw <- map(pop_nw, ~sim_survey(.$pop,
#                                    n_sims = nsurveys, # one survey per item in population list object
#                                    trawl_dim = trawl_dim,
#                                    q = catch_q,
#                                    set_den = set_den,
#                                    min_sets = min_sets,
#                                    age_sampling = age_sampling,
#                                    age_space_group = age_space_group,
#                                    resample_cells = resample_cells))
#

## SAVE THE DATA ####
saveRDS(survdat_sq, here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat.rds", sep = "_")))
#saveRDS(survdat_sq_nw, here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat_nw.rds", sep = "_")))


