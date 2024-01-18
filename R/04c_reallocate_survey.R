### created: 01/18/2024
### updated:

# 04c - SIMULATE REALLOCATED SURVEY ####


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


### DATA SET UP ####
# data locations
sdmtmb.dir <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis/sdmtmb"
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

### number of simulations
nsims <- 1:2

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

# simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
survdat_sq <- readRDS(here(survdat, str_c(species, season, "sq-survdat.rds", sep = "_")))

# simulated status quo survey data created here("R", "04b_preclude_survey.R")
survdat_precl <- readRDS(here(survdat, str_c(species, season, "precl-survdat.rds", sep = "_")))


## Identify tows ####
#tows occuring inside wind areas | AREA CODE = 1
wind_tows <- map(survdat_sq, ~filter(.$setdet, AREA_CODE == 1))

wind_cells <- map(wind_tows, ~group_by(.,sim,year,strat) |>
                    distinct(cell)
)

wind_strat <- map(wind_cells, ~unique(.$strat))


## Filter grid ####
# extract grids
grids <- map(survdat_sq, ~pluck(., "grid_xy"))

# filter grids for strata that have wind but cells that are outside wind areas
out_wa_grid <- map2(grids, wind_strat, ~filter(.x, strat %in% .y, AREA_CODE == 2) |> group_by(strat) |> nest())

#extract unique strata that have wind overlaps
out_strat <- map(out_wa_grid, ~unique(.$strat))

# count the number of wind tows in each year and strata that need reallocating
wind_summ <- map(wind_tows, ~group_by(.x, sim, year, strat) |>
                   nest() |>
                   mutate(count = map(data, ~length(.$set))) |> rename(wind_tows = data))

## Generate new locations ####
# join the count data with the grid in order to sample the grip
join_data <- map2(wind_summ, out_wa_grid, ~left_join(.x, .y, by="strat"))

#filter out the nulls (strata where no cells remain outside of wind areas)
join_data2 <- map2(join_data, out_strat, ~filter(.x, strat %in% .y))
# check that only nulls were removed
map2(join_data, join_data2, ~anti_join(.x, .y, by=c("sim", "strat", "year")))

# sample new locations and create a set dataframe to supply to custom_sets argument
new_locations <- map(join_data2,
                     ~mutate(., new = map2(data, count,
                                                ~slice_sample(.x, n=.y, replace=TRUE)),
                             tow_info = map(wind_tows, ~select(., !c(x,y,cell,depth, AREA_CODE, n, n_aged, n_measured, N))),
                             new_set_loc = map2(new, tow_info, ~bind_cols(.x, .y))) |>
                       select(c(sim, year, strat, new_set_loc)) |>
                       unnest(cols=new_set_loc) |>
                       #select(!(division)) |>
                       relocate(set, .after = last_col()) |>
                       select(!c(AREA, division)))

# less tows bc could not reallocate 0320 tows

## Simulate New Tow Data ####
set.seed(13487)

survdat_new_locs <- map2(pop, new_locations, ~sim_survey(.x$pop,
                                   n_sims = 1, # one survey per item in population list object
                                   trawl_dim = trawl_dim,
                                   q = catch_q,
                                   set_den = set_den,
                                   min_sets = min_sets,
                                   age_sampling = age_sampling,
                                   age_space_group = age_space_group,
                                   custom_sets = .y,
                                   resample_cells = resample_cells)
)

## Bind to precluded survey ####
survdat_reall <- map2(survdat_precl, survdat_new_locs, ~bind_rows(.x, .y$setdet))


### to do: add code to bind survdat_new_locs$samps for length and age data to precluded survey for survey data product comps calculations


## SAVE THE DATA ####
saveRDS(survdat_new_locs, here(survdat, str_c(species, season, "new-locs-survdat.rds", sep = "_")))
saveRDS(survdat_reall, here(survdat, str_c(species, season, "reall-survdat.rds", sep = "_")))


