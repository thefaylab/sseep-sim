### created: 01/18/2024
### updated: 02/07/2025

# 04c - SIMULATE REALLOCATED SURVEY ####


## Objective ####
# identify and calculate the number of sets that occurred in a wind area and the strata in which those sets and wind areas occurred
# remove wind cells from the grid and sample the remaining cells in the same strata that precluded wind tows occurred and for the same number of tows that were precluded
# set these new sampled cells as new locations and simulate a status quo NMFS bottom trawl survey.
# bind the resulting tow level data to the dataset where wind tows and catch rates were removed

# Outputs: locations of reallocated tows in strata where wind tows were precluded, and a full survey with tow level data for each replicate of a simulated population and abundance.

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
library(tidyverse)
library(tidyr)
library(purrr)
library(dplyr)
# source(here("R", "sim_pop_fn.R"))
set.seed(131)

### DATA SET UP ####
# data locations
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
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
source(("R/selectivity_fns.R"))

#catch_q <- sim_logistic(k = 2, x0 = 2.5) #older sel
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
#pop <- readRDS(here(dist.dat, str_c(species, season, "50abund-dist.rds", sep = "_")))
pop <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "abund-dist.rds", sep = "_")))

# pop_nw <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "nowind-abund-dist.rds", sep = "_"))) # with wind effect turned off

# simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
survdat_sq <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat.rds", sep = "_")))
# survdat_sq_nw <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat_nw.rds", sep = "_"))) # where wind effect is turned off

# simulated status quo survey data created here("R", "04b_preclude_survey.R")
survdat_precl <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "precl-surv-dat.rds", sep = "_")))
# survdat_precl_nw <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "precl-surv-dat_nw.rds", sep = "_")))
 # where wind effect is turned off

## Identify tows ####
#tows occuring inside wind areas | AREA CODE = 1
#wind_tows <- map(survdat_sq, ~filter(.$setdet, AREA_CODE == 1))
wind_tows <- map(survdat_sq, ~filter(.$setdet, AREA_CODE == 1, year >= 6))


# wind_cells <- map(wind_tows, ~group_by(.,sim,year,strat) |>
#                     distinct(cell))

wind_cells <- map(wind_tows, ~group_by(., sim, year, strat) |>
                    filter(year >= 6) |> # Apply change only for years 6-15
                    distinct(cell))


wind_strat <- map(wind_cells, ~unique(.$strat))


## Filter grid ####
# extract grids
grids <- map(survdat_sq, ~pluck(., "grid_xy"))

# filter grids for strata that have wind but cells that are outside wind areas
out_wa_grid <- map2(grids, wind_strat, ~filter(.x, strat %in% .y, AREA_CODE == 2) |> group_by(strat) |> nest())

#extract unique strata that have wind overlaps
out_strat <- map(out_wa_grid, ~unique(.$strat))

# count the number of wind tows in each year and strata that need reallocating
# wind_summ <- map(wind_tows, ~group_by(.x, sim, year, strat) |>
#                    nest() |>
#                    mutate(count = map(data, ~length(.$set))) |> rename(wind_tows = data))

wind_summ <- map(wind_tows, ~group_by(.x, sim, year, strat) |>
                   filter(year >= 6) |> # Apply only for years 6-15
                   nest() |>
                   mutate(count = map(data, ~length(.$set))) |> 
                   rename(wind_tows = data))


## Generate new locations ####
# join the count data with the grid with remaining cells outside wind areas in order to sample the grid
join_data <- map2(wind_summ, out_wa_grid, ~left_join(.x, .y, by="strat"))

#filter out the nulls (strata where no cells remain outside of wind areas)
join_data2 <- map2(join_data, out_strat, ~filter(.x, strat %in% .y))

# check that only nulls were removed
map2(join_data, join_data2, ~anti_join(.x, .y, by=c("sim", "strat", "year"))) |> head(2)

# sample new locations and create a set dataframe to supply to custom_sets argument
# new_locations <- map(join_data2,
#                      ~mutate(., new = purrr::map2(data, count,
#                                                 ~slice_sample(.x, n=.y, replace=TRUE)),
#                              tow_info = purrr::map(wind_tows,
#                                                    ~select(., !c(x,y,cell,depth,AREA_CODE, n, n_aged, n_measured, N))),
#                              new_set_loc = purrr::map2(new, tow_info, ~bind_cols(.x, .y))) |>
#                        dplyr::select(c(sim, year, strat, new_set_loc)) |>
#                        tidyr::unnest(cols=new_set_loc) |>
#                        #select(!(division)) |>
#                        dplyr::relocate(set, .after = last_col()) |>
#                        dplyr::select(!c(AREA)))


new_locations <- map(join_data2,
                     ~mutate(., new = purrr::map2(data, count,
                                                  ~slice_sample(.x |> filter(year >= 6), n=.y, replace=TRUE)), # Apply change only for years 6-15
                             tow_info = purrr::map(wind_tows,
                                                   ~select(., !c(x, y, cell, depth, AREA_CODE, n, n_aged, n_measured, N))),
                             new_set_loc = purrr::map2(new, tow_info, ~bind_cols(.x, .y))) |>
                       dplyr::select(c(sim, year, strat, new_set_loc)) |>
                       tidyr::unnest(cols=new_set_loc) |>
                       dplyr::relocate(set, .after = last_col())) #|>
                     #  dplyr::select(!c(AREA)))


#Check
map(new_locations, ~unique(.x$year)) #should contain only years 6 to 15
map(new_locations, ~count(.x, year)) #ensure that only years 6-15 contain new data. If years 1-5 are missing from the output, it confirms the fix
map(new_locations, ~filter(.x, year < 6)) #To double-check that years 1-5 remain unchanged, an empty list confirms that no changes were made to years 1-5
#To confirm that only years 6-15 changed, compare the old and new tow allocations
map2(join_data2, new_locations, ~full_join(count(.x, year), count(.y, year), by = "year", suffix = c("_before", "_after")))


# less tows bc could not reallocate 0320 tows

## Simulate New Tow Data ####
survdat_new_locs <- map2(pop, new_locations, ~sim_survey(.x$pop,
                                   n_sims = nsurveys, # one survey per item in population list object
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


map(survdat_reall, ~unique(.x$year))
map(survdat_reall, ~table(.x$year, .x$AREA_CODE))
map(survdat_reall, ~table(.x$sim, .x$AREA_CODE, .x$year))

### to do: add code to bind survdat_new_locs$samps for length and age data to precluded survey for survey data product comps calculations



library(ggplot2)

map(survdat_reall, ~ggplot(.x, aes(x = as.factor(year), fill = as.factor(AREA_CODE))) +
      geom_bar() +
      theme_minimal() +
      labs(title = "Survey Tow Distribution After Reallocation", 
           x = "Year", y = "Number of Tows", fill = "AREA_CODE"))


## SAVE THE DATA ####
saveRDS(survdat_new_locs, here(survdat, str_c(species, season, length(nsims), "sims", "locs", nsurveys, "survdat.rds", sep = "_")))
saveRDS(survdat_reall, here(survdat, str_c(species, season, length(nsims), "sims", "reall", nsurveys, "survdat.rds", sep = "_")))
# saveRDS(survdat_new_locs, here(survdat, str_c(species, season, length(nsims), "sims", "locs", nsurveys, "survdat_nw.rds", sep = "_")))
# saveRDS(survdat_reall, here(survdat, str_c(species, season, length(nsims), "sims", "reall", nsurveys, "survdat_nw.rds", sep = "_")))
