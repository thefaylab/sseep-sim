### created: 01/18/2024
### updated: 04/16/2025

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
library(tidyverse)
library(data.table)
library(here)
source(here("R", "selectivity_fns.R"))
set.seed(123)

### DATA SET UP ####
# Directories
dist.dat  <- here("data", "rds", "dists")
survdat   <- here("data", "rds", "survdat")

# Parameters
species   <- "scup"
season    <- "fall"
ages      <- 0:7
years     <- 1:15
nsims     <- 1:100
nsurveys  <- 25
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))


# Trawl survey configuration
trawl_dim       <- c(2.7, 0.014)
catch_q         <- force_sim_logistic(k = -0.66, x0 = -1.14, force_age = TRUE, age = 0, force_sel = 1)
set_den         <- 0.001
min_sets        <- 3
age_sampling    <- "stratified"
age_space_group <- "set"
resample_cells  <- TRUE


# Set which chunk to run (1 to 5)
chunk_id <- 1  # Change this (1 to 5)
this_chunk <- chunks[[chunk_id]]


# LOOP 1: Generate new tow locations
for (i in this_chunk) {
  message(sprintf("Generating reallocated locations for population %03d of chunk %d...", i, chunk_id))

  # Load SQ survey
  sq_survey <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_sq_survey.rds", species, season, i, nsurveys)))
  grid <- sq_survey$grid_xy

  # Identify tows in wind areas after year 5
  wind_tows <- sq_survey$setdet |> dplyr::filter(AREA_CODE == 1, year >= 6)
  wind_cells <- wind_tows |> dplyr::group_by(sim, year, strat) |> dplyr::distinct(cell)
  wind_strat <- unique(wind_cells$strat)

  # Filter grid outside wind areas for affected strata
  grid_outside <- grid |> dplyr::filter(AREA_CODE == 2, strat %in% wind_strat)
  out_wa_grid <- grid_outside |> dplyr::group_by(strat) |> tidyr::nest()
  out_strat <- unique(out_wa_grid$strat)

  # Count wind tows to reallocate
  wind_summ <- wind_tows |>
    dplyr::group_by(sim, year, strat) |>
    tidyr::nest() |>
    dplyr::mutate(count = purrr::map(data, ~length(.$set))) |>
    dplyr::rename(wind_tows = data)

  # Join with grid and remove strata with no remaining area
  join_data <- dplyr::left_join(wind_summ, out_wa_grid, by = "strat") |>
    dplyr::filter(strat %in% out_strat) |>
    dplyr::filter(!purrr::map_lgl(data, is.null))

  # Generate new locations
  new_locations <- join_data |>
    dplyr::mutate(
      new = purrr::map2(data, count, ~dplyr::slice_sample(.x, n = .y, replace = TRUE)),
      new_set_loc = purrr::pmap(list(new, wind_tows, sim), ~{
        tow_info <- ..2 |> dplyr::filter(sim == ..3) |>
          dplyr::select(-c(x, y, cell, depth, AREA_CODE, n, n_aged, n_measured, N))
        dplyr::bind_cols(..1, tow_info)
      })
    ) |>
    dplyr::select(sim, year, strat, new_set_loc) |>
    tidyr::unnest(cols = new_set_loc) |>
    dplyr::relocate(set, .after = dplyr::last_col())

  # Save output
  saveRDS(new_locations, here(survdat, sprintf("%s_%s_%03d_%d_locs_survey.rds", species, season, i, nsurveys)))
}




 #Check
new_locs_084 <- readRDS(here(survdat, "scup_fall_084_25_locs_survey.rds"))
old_locs <- readRDS(here(survdat, "scup_fall_2_sims_locs_25_survdat.rds"))
old_locs[[1]]$setdet
(new_locs_084)

# Tows per year
new_locs_084 |> count(year)

# Tows per year and stratum
new_locs_001 |> count(year, strat)

# Visualize them
library(ggplot2)
ggplot(new_locs_001, aes(x, y)) +
  geom_point(aes(color = as.factor(year)), alpha = 0.6) +
  facet_wrap(~year) +
  theme_minimal() +
  labs(title = "Reallocated Tow Locations for Pop 001")



#LOOP 2: Generate new reallocated survey
set.seed(123)
chunk_id <- 5  # Change this (1 to 5)
this_chunk <- chunks[[chunk_id]]

for (i in this_chunk) {
  message(sprintf("Running reallocated survey for population %03d of chunk %d...", i, chunk_id))

  # Load population and precluded survey
  pop_i         <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))
  precl_survey  <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_precl_survey.rds", species, season, i, nsurveys)))
  new_locations <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_locs_survey.rds", species, season, i, nsurveys)))

  # Run survey with reallocated sets
  survey_new <- sim_survey(
    pop_i,
    n_sims           = nsurveys,
    trawl_dim        = trawl_dim,
    q                = catch_q,
    set_den          = set_den,
    min_sets         = min_sets,
    age_sampling     = age_sampling,
    age_space_group  = age_space_group,
    custom_sets      = new_locations,
    resample_cells   = resample_cells
  )

  # Combine with precluded survey
  survey_reall <- bind_rows(precl_survey, survey_new$setdet)

  # Save outputs
  saveRDS(new_locations, here(survdat, sprintf("%s_%s_%03d_%d_locs_survey.rds", species, season, i, nsurveys)))
  saveRDS(survey_reall, here(survdat, sprintf("%s_%s_%03d_%d_reall_survey.rds", species, season, i, nsurveys)))
}


survdat_reall_old <- readRDS(here(survdat, "scup_fall_2_sims_reall_25_survdat.rds"))
survdat_reall_1 <- readRDS(here(survdat, "scup_fall_001_25_reall_survey.rds"))
survdat_reall_1




survdat_reall_old[[1]]




#Count
s_real <- survdat_reall_1

# Count sets by sim, year
tow_counts3 <- s_real |>
  count(sim, year, name = "n_sets") #can add strat to check

print(tow_counts3)


#Average number of tows per year across all simulations
avg_tows_per_year3 <- s_real |>
  count(sim, year) |>
  group_by(year) |>
  summarise(
    mean_tows = mean(n),
    sims = n(),
    .groups = "drop"
  )

print(avg_tows_per_year3)




#Check which strata are affected
strata_with_wind <- unique(wind_summ$strat)
strata_with_available_cells <- unique(out_wa_grid$strat)
strata_failed <- setdiff(strata_with_wind, strata_with_available_cells)
strata_failed



# ### LOAD DATA ####
# # simulated abundance and distributions
# pop <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "abund-dist.rds", sep = "_")))
#
# # simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
# survdat_sq <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat.rds", sep = "_")))
#
# # simulated status quo survey data created here("R", "04b_preclude_survey.R")
# survdat_precl <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "precl-surv-dat.rds", sep = "_")))
#
# ## Identify tows ####
# #tows occuring inside wind areas | AREA CODE = 1
# wind_tows <- map(survdat_sq, ~filter(.$setdet, AREA_CODE == 1, year >= 6))
#
# wind_cells <- map(wind_tows, ~group_by(., sim, year, strat) |>
#                     filter(year >= 6) |> # Apply change only for years 6-15
#                     distinct(cell))
#
# wind_strat <- map(wind_cells, ~unique(.$strat))
#
#
# ## Filter grid ####
# # extract grids
# grids <- map(survdat_sq, ~pluck(., "grid_xy"))
#
# # filter grids for strata that have wind but cells that are outside wind areas
# out_wa_grid <- map2(grids, wind_strat, ~filter(.x, strat %in% .y, AREA_CODE == 2) |> group_by(strat) |> nest())
#
# #extract unique strata that have wind overlaps
# out_strat <- map(out_wa_grid, ~unique(.$strat))
#
# # count the number of wind tows in each year and strata that need reallocating
# wind_summ <- map(wind_tows, ~group_by(.x, sim, year, strat) |>
#                    filter(year >= 6) |> # Apply only for years 6-15
#                    nest() |>
#                    mutate(count = map(data, ~length(.$set))) |>
#                    rename(wind_tows = data))
#
#
# ## Generate new locations ####
# # join the count data with the grid with remaining cells outside wind areas in order to sample the grid
# join_data <- map2(wind_summ, out_wa_grid, ~left_join(.x, .y, by="strat"))
#
# #filter out the nulls (strata where no cells remain outside of wind areas)
# join_data2 <- map2(join_data, out_strat, ~filter(.x, strat %in% .y))
#
# # check that only nulls were removed
# map2(join_data, join_data2, ~anti_join(.x, .y, by=c("sim", "strat", "year"))) |> head(2)
#
#
# new_locations <- map(join_data2,
#                      ~mutate(., new = purrr::map2(data, count,
#                                                   ~slice_sample(.x |> filter(year >= 6), n=.y, replace=TRUE)), # Apply change only for years 6-15
#                              tow_info = purrr::map(wind_tows,
#                                                    ~select(., !c(x, y, cell, depth, AREA_CODE, n, n_aged, n_measured, N))),
#                              new_set_loc = purrr::map2(new, tow_info, ~bind_cols(.x, .y))) |>
#                        dplyr::select(c(sim, year, strat, new_set_loc)) |>
#                        tidyr::unnest(cols=new_set_loc) |>
#                        dplyr::relocate(set, .after = last_col()))
#
#
# #Check
# map(new_locations, ~unique(.x$year)) #should contain only years 6 to 15
# map(new_locations, ~count(.x, year)) #ensure that only years 6-15 contain new data. If years 1-5 are missing from the output, it confirms the fix
# map(new_locations, ~filter(.x, year < 6)) #To double-check that years 1-5 remain unchanged, an empty list confirms that no changes were made to years 1-5
# #To confirm that only years 6-15 changed, compare the old and new tow allocations
# map2(join_data2, new_locations, ~full_join(count(.x, year), count(.y, year), by = "year", suffix = c("_before", "_after")))
#
#
# # less tows bc could not reallocate 0320 tows
#
# ## Simulate New Tow Data ####
# survdat_new_locs <- map2(pop, new_locations, ~sim_survey(.x$pop,
#                                    n_sims = nsurveys, # one survey per item in population list object
#                                    trawl_dim = trawl_dim,
#                                    q = catch_q,
#                                    set_den = set_den,
#                                    min_sets = min_sets,
#                                    age_sampling = age_sampling,
#                                    age_space_group = age_space_group,
#                                    custom_sets = .y,
#                                    resample_cells = resample_cells)
# )
#
# ## Bind to precluded survey ####
# survdat_reall <- map2(survdat_precl, survdat_new_locs, ~bind_rows(.x, .y$setdet))
#
#
# map(survdat_reall, ~unique(.x$year))
# map(survdat_reall, ~table(.x$year, .x$AREA_CODE))
# map(survdat_reall, ~table(.x$sim, .x$AREA_CODE, .x$year))
#
# ### to do: add code to bind survdat_new_locs$samps for length and age data to precluded survey for survey data product comps calculations
#
#
#
# library(ggplot2)
#
# map(survdat_reall, ~ggplot(.x, aes(x = as.factor(year), fill = as.factor(AREA_CODE))) +
#       geom_bar() +
#       theme_minimal() +
#       labs(title = "Survey Tow Distribution After Reallocation",
#            x = "Year", y = "Number of Tows", fill = "AREA_CODE"))
#
#
# ## SAVE THE DATA ####
# saveRDS(survdat_new_locs, here(survdat, str_c(species, season, length(nsims), "sims", "locs", nsurveys, "survdat.rds", sep = "_")))
# saveRDS(survdat_reall, here(survdat, str_c(species, season, length(nsims), "sims", "reall", nsurveys, "survdat.rds", sep = "_")))
# # saveRDS(survdat_new_locs, here(survdat, str_c(species, season, length(nsims), "sims", "locs", nsurveys, "survdat_nw.rds", sep = "_")))
# # saveRDS(survdat_reall, here(survdat, str_c(species, season, length(nsims), "sims", "reall", nsurveys, "survdat_nw.rds", sep = "_")))
