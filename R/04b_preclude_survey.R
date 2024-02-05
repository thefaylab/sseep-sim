### created: 01/18/2024
### updated: 02/05/2024

# 04b - PRECLUDE SURVEY FROM WIND ENERGY AREAS ####


## Objective ####
# For a given species, distribution, and survey, sets are removed from the simulated tow level data if they occurred in a cell indexed as a wind cell representing an overlap with offshore wind energy areas.

# Outputs: one survey and respective tow level data occurring outside of offshore wind energy areas only for each replicate of a simulated population and abundance

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
# source(here("R", "sim_pop_fn.R"))


### DATA SET UP ####
# data locations
sdmtmb.dir <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis/sdmtmb"
survdat <- here("data", "rds", "survdat")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"



### LOAD DATA ####
# simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
survdat_sq <- readRDS(here(survdat, str_c(species, season, "sq-survdat.rds", sep = "_")))


## PRECLUDE SURVEY ####
# remove tows/sets that were sampled in cells indexed as a wind cell to represent the preclusion of the survey due to offshore wind energy areas
# `AREA_CODE == 2` represents outside areas
survdat_precl <- map(survdat_sq, ~filter(.$setdet, AREA_CODE == 2))


### to do: add code to filter survdat_sq$samps for length and age data calculations for samples only taken outside of wind energy areas


## SAVE THE DATA ####
saveRDS(survdat_precl, here(survdat, str_c(species, season, "precl-surv-dat.rds", sep = "_")))
