### created: 01/18/2024
### updated:

# 03 - APPEND DISTRIBUTIONS ####


## Objective ####
# For a given species, distribute the simulated abundances generated here("R", "02_simulate_populations.R") across the bottom trawl survey footprint using the spatial distribution predictions from a spatiotemporal generalized linear mixed model fit to historic survey observations.

# Outputs: simulated abundance and distribution for across a set of iterations

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
pop.dat <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")
dist.dat <- here("data", "rds", "dists")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### number of simulations
nsims <- 1:50

### LOAD DATA ####
### the NE bts survey grid
# as dataframe
grid_xy <- readRDS(here("data", "rds", "survey_grid_062022.rds")) |>
  rename(strat = STRATUM,
         x = X,
         y = Y) |>#, # rename to match column calls within SimSurvey functions
  #depth = AVGDEPTH) |>
  mutate(division = 1) |>#, # add division information
  #x = X/1000, # convert to km
  #y = Y/1000) |>
  # dplyr::select(x, y, cell, division, strat, depth, AREA_CODE) |>
  data.table::as.data.table()

# as stars object
grid_stars <- readRDS(here("data",  "survey_grid_stars_062022.rds"))

### species seasonal spatial distribution predictions
preds <- readRDS(file = here(sdmtmb.dir, species, "data", str_c(season, "_grid_preds.rds", sep = "")))

# load spring model predictions
# spring_preds <- readRDS(file = here(sdmtmb.dir, species, "data", "spring_grid_preds.rds"))

### simulated species abundance
pop <- readRDS(here(pop.dat, str_c(species, "_pop.rds", sep = "")))

### simulated species N at age matrices
Nage <- readRDS(here(Nage.dat, str_c(species, "_Nage.rds", sep = "")))


## Distribute numbers at age ####
### Filter predictions ####
preds <- preds |>
  filter(EST_YEAR %in% c(2014:2016, 2018, 2019)) |>
  rename(year = EST_YEAR, #rename to match column calls in SimSurvey
         strat = STRATUM,
         depth = AVGDEPTH,
         x = X,
         y = Y) |>
  mutate(N_dist = exp(est),
         year = case_when( # change the year values based on their sequence
           year == 2014 ~ 1,
           year == 2015 ~ 2,
           year == 2016 ~ 3,
           year == 2018 ~ 4,
           year == 2019 ~ 5
         )) |>#,
  #cell = seq(1:length(N_dist))) |> # add cell # value
  dplyr::select(x,y, year, N_dist, cell, strat, depth) |>
  data.table::as.data.table()

### Calculate probability of distributions ####
dist <- sdmTMB::replicate_df(preds, "age", ages) |> # replicate the predictions over each age, so there are distributions for each age
  sdmTMB::replicate_df("sim", nsims) |> # replicate the predictions over the number of iterations, so there are distributions for each replicate
  left_join(Nage, by = c("sim", "age")) |> # join populations at age to the predictions
  rename(Nage = "pluck_n0") |>
  group_by(year, age, sim) |>
  mutate(P_i = N_dist/sum(N_dist), # probability of distribution
         N = Nage * P_i, # multiply probability by simulated numbers of age
         age = as.double(age),
         cell = as.double(cell),
         division = 1) |>
  dplyr::select(sim, age, year, cell, N, x, y, strat, division, depth)

### Force zeros at age ####
no_dist <- dplyr::anti_join(grid_xy, preds, by = c("x", "y", "cell")) |>
  sdmTMB::replicate_df("year", years) |>
  sdmTMB::replicate_df("age", ages) |>
  sdmTMB::replicate_df("sim", nsims) |>
  mutate(N = 0,
         age = as.double(age),
         cell = as.double(cell),
         division = 1)

### Bind distributions ####
full_dist <- bind_rows(no_dist, dist)

## Append population object ####
### The grid ####
# create list object with one item containing both grid objects
grids <- list(list(grid = grid_stars, grid_xy = grid_xy))

# replicate the grid list over the same number of iterations so each list item replicate contains both grids
grids <- rep(grids, length(nsims))

# append the replicated grids list to each iterated simulated abundance
pop <- map2(pop, grids, ~append(.x, .y))

### Spatial distribution ####
# create list object from full distribution object with each item derived by replicate
dist_list <- full_dist |>
  split(by = "sim") |>
  map(~list(sp_N = .)) # title each item as `sp_N` to match SimSurvey calls

# append the distributions to the simulated abundance object
pop <- map2(pop, dist_list, ~append(.x, .y)) |>
  map(~list(pop = .)) # prep for mapping through sim_survey() in future script

## SAVE THE DATA ####
saveRDS(pop, here(dist.dat, str_c(species, season, "50abund-dist.rds", sep = "_")))


