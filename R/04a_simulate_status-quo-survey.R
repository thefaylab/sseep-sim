### created: 01/18/2024
### updated: 02/05/2024

# 04a - SIMULATE STATUS QUO SURVEY ####


## Objective ####
# For a given species and distribution, simulate the status quo NMFS bottom trawl survey.

# Outputs: one survey and respective tow level data for each replicate of a simulated population and abundance

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
# source(here("R", "sim_pop_fn.R"))
set.seed(1462)


### DATA SET UP ####
# data locations
sdmtmb.dir <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis/sdmtmb"
pop.dat <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### trawl dimensions
trawl_dim <- c(2.7, 0.014)

### catchability from the most recent stock assessment
catch_q <- sim_logistic(k = 2, x0 = 2.5)

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
pop <- readRDS(here(dist.dat, str_c(species, "_abund-dist.rds", sep = "")))



## SIMULATE SURVEY ####
set.seed(13487)

survdat_sq <- map(pop, ~sim_survey(.$pop,
                                   n_sims = 1, # one survey per item in population list object
                                   trawl_dim = trawl_dim,
                                   q = catch_q,
                                   set_den = set_den,
                                   min_sets = min_sets,
                                   age_sampling = age_sampling,
                                   age_space_group = age_space_group,
                                   resample_cells = resample_cells)
)


## SAVE THE DATA ####
saveRDS(survdat_sq, here(survdat, str_c(species, season, "sq-surv-dat.rds", sep = "_")))


