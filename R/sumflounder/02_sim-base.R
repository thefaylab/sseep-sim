## created: 11/03/2023
### updated:

# 02 - SIMULATE SS-BASE SCENARIO ####

## Objective ####
# Script will:
## simulate 1000 summer flounder population and abundances based on stock assessment data
## generate 1000 summer flounder distributions based on sdmTMB predictions and simulated populations
## simulate 1000 surveys on each population
##
#
#

### LOAD PACKAGES ####
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(patchwork)
library(here)
library(janitor)
library(sdmTMB)
set.seed(380)
theme_set(theme_bw())


sseep.analysis <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis"
sdmtmb.dir <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis/sdmtmb"

### DATA ####
grid_xy <- readRDS(here("data", "rds", "survey_grid_062022.rds")) |>
  rename(strat = STRATUM,
         x = X,
         y = Y) |>#, # rename to match column calls within SimSurvey functions
  #depth = AVGDEPTH) |>
  mutate(division = 1) |>#, # add division information
  #x = X/1000, # convert to km
  #y = Y/1000) |>
  dplyr::select(x, y, cell, division, strat, depth) |>
  data.table::as.data.table()

# the survey Grid as stars object
grid_stars <- readRDS(here("data",  "survey_grid_stars_062022.rds"))

# load spring model predictions
spring_preds <- readRDS(file = here(sdmtmb.dir, "sumflounder", "data", "spring_predictions.rds"))

# Recruitment at age 0 for the most recent year, 2016 and 2017; numbers were reported in thousands in the 66th SAW report.
Rec_age0 <- c(43000, 44552)*1000

# fishing mortality at age for the most recent year, 2016 and 2017
ages <- as.character(0:7)
years <- as.character(1:2)
F <- matrix(c(0.011, 0.045, 0.127, 0.253, 0.417, 0.388, 0.381, 0.277, 0.009, 0.043, 0.115, 0.213, 0.334, 0.303, 0.295, 0.217), nrow = 8, ncol = 2, byrow = FALSE, dimnames = list(age = 0:7, year = 1:2))

# natural mortality
M <- 0.25

#total mortality
Z <- F + M

# Von Bertalanffy growth parameters for both male and female
Linf = 83.6
K = 0.14


## SIMULATE ABUNDANCE ####
# simulate a summer flounder like population 10 times
pop <- map(1:10, ~sim_abundance(ages = 0:7, years = 1:2,
                     R = sim_R(log_mean = log(Rec_age0), log_sd = 0.001, plot = TRUE),
                     Z = sim_Z(log_mean = log(Z), log_sd = 0.001, plot = TRUE),
                     N0 = sim_N0(N0 = "exp", plot = TRUE),
                     growth = sim_vonB(Linf = 83.6, K = 0.14, plot = TRUE)))

# save each initial numbers at age (N0) object from each simulated population
Nage <- map_df(pop, ~pluck(.,"N0")|> tibble(age = 0:7), .id = "replicate") |>
  clean_names() |>
  mutate(replicate = as.integer(replicate)) |>
  rename(N0 = pluck_n0)

# save the data
saveRDS(pop, here("data", "trials", "pop-sim10.rds"))
saveRDS(Nage, here("data", "trials", "nage-sim10.rds"))

## APPEND GRIDS ####
# create a grids list object that contains both the grid data table and starts object to append to the simulated population object
grids <- list(list(grid_xy = grid_xy, grid = grid_stars))

# replicate list object 10 times to match the number of population replicates
grids <- rep(grids, 10)

# append each population replicate with the grid object
pop <- map2(pop, grids, ~append(.x, .y))

## APPREND DISTRIBUTION ####
# extract the predictions from sdmTMB
preds_2yr <- spring_preds |>
  filter(EST_YEAR %in% c(2016,2017)) |>
  rename(year = EST_YEAR, #rename to match column calls in SimSurvey
         strat = STRATUM,
         depth = AVGDEPTH) |>
  mutate(N_dist = exp(est),
         year = case_when( # change the year values based on their sequence
           year == 2016 ~ 1,
           year == 2017 ~ 2
         )) |>#,
  #cell = seq(1:length(N_dist))) |> # add cell # value
  dplyr::select(X,Y, year, N_dist, cell, strat, depth) |>
  data.table::as.data.table()

# force 0s at age outside sdmTMB prediction area
no_dist <- dplyr::anti_join(grid_xy, preds_2yr, by = "cell") |>
  sdmTMB::replicate_df("age", 0:7) |>
  sdmTMB::replicate_df("year", 1:2) |>
  sdmTMB::replicate_df("replicate", 1:10) |>
  mutate(N = 0, #
         age = as.double(age),
         cell = as.double(cell),
         division = 1)

# predicted numbers at age based on simulated N0 and sdmtmb predictions
dist <- sdmTMB::replicate_df(preds_2yr, "age", 0:7) |> # replicate the predictions over each age, so there are distributions for each age
  sdmTMB::replicate_df("replicate", 1:10) |>
  left_join(Nage, by = c("replicate", "age")) |> # join populations at age to the predictions
  rename(x = X,
         y = Y) |>
  group_by(replicate) |>
  mutate(P_i = N_dist/sum(N_dist), # probability of distribution
         N = N0 * P_i, # multiply probability by simulated numbers of age
         age = as.double(age),
         cell = as.double(cell),
         division = 1) |>
  dplyr::select(age, year, cell, N, x, y, strat, replicate, depth, division)

# bind the distributed numbers at age to cover the full spatial extent of NEUS
full_dist <- bind_rows(dist, no_dist)

# pull each replicate into a list object
dist_list <- full_dist |>
  split(~as.factor(replicate))

# name each list object sp_N to match data calls in SimSurvey
dist_list <- map(dist_list, ~list(sp_N = .))

# append each distribution object to the respective simulated population
pop <- map2(pop, dist_list, ~append(.x, .y))

# name each simulated pop and distribution object pop to match the data calls in SimSurvey
pop <- map(pop, ~list(pop = .))

# save the data
saveRDS(pop, here("data", "trials", "dist-sim10.rds"))

## SIMULATE SURVEY ####
# simulate a single survey over each of the ten populations
surv_dat <- map(pop, ~sim_survey(.$pop, n_sims = 1,
                                 trawl_dim = c(2.7, 0.014),
                                 q = sim_logistic(k = 2, x0 = 2.5), #if plot = TRUE, xlim errors
                                 set_den = 0.001,
                                 min_sets = 3,
                                 #lengths_cap = # max # of lengths to collect per set
                                 #age_cap = 8, # max # of ages to sample per length group - all 8?
                                 age_sampling = "stratified",
                                 #age_length_group = , # length group bin size (cm) -
                                 age_space_group = "set", # spatial scale of stratified age sampling - ages taken each tow?
                                 resample_cells = TRUE))

# save the data
saveRDS(surv_dat, here("data", "trials", "survdat-10sims.rds"))
