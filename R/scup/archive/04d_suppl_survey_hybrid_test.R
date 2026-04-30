### created: 10/14/2025
### updated:

# 04c - SIMULATE SUPPLEMENTAL SAMPLING ####
# random stratified outside
# fixed stations inside wind areas


## Objective ####

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
age_sampling    <- "random"
age_space_group <- "set"
resample_cells  <- TRUE


# Set which chunk to run (1 to 5)
chunk_id <- 1  # Change this (1 to 5)
this_chunk <- chunks[[chunk_id]]


### LOOP 1: Generate random tow locations INSIDE wind areas (fixed across years 6–15)
for (i in this_chunk) {
  message(sprintf("Generating inside-wind random locations for population %03d of chunk %d...", i, chunk_id))

  # Load SQ survey
  sq_survey <- readRDS(here(
    survdat,
    sprintf("%s_%s_%03d_%d_sq_survey.rds", species, season, i, nsurveys)
  ))
  grid <- sq_survey$grid_xy

  # Identify tows inside wind areas after year 5
  wind_tows <- sq_survey$setdet |>
    dplyr::filter(AREA_CODE == 1, year >= 6)
  wind_cells <- wind_tows |>
    dplyr::group_by(sim, year, strat) |>
    dplyr::distinct(cell)
  wind_strat <- unique(wind_cells$strat)

  # Grid restricted to INSIDE wind areas
  grid_inside <- grid |>
    dplyr::filter(AREA_CODE == 1, strat %in% wind_strat)
  in_wa_grid <- grid_inside |>
    dplyr::group_by(strat) |>
    tidyr::nest()
  in_strat <- unique(in_wa_grid$strat)

  # Count number of tows to reallocate (year 6 only)
  wind_summ <- wind_tows |>
    dplyr::filter(year == 6) |>
    dplyr::group_by(sim, year, strat) |>
    tidyr::nest() |>
    dplyr::mutate(count = purrr::map(data, ~length(.$set))) |>
    dplyr::rename(wind_tows = data)

  # Join with grid and remove strata with no available cells
  join_data <- dplyr::left_join(wind_summ, in_wa_grid, by = "strat") |>
    dplyr::filter(strat %in% in_strat) |>
    dplyr::filter(!purrr::map_lgl(data, is.null))

  # Generate new locations (sample once inside wind areas)
  new_locations_y6 <- join_data |>
    dplyr::mutate(
      new = purrr::map2(data, count,
                        ~dplyr::slice_sample(.x, n = .y, replace = TRUE)),
      new_set_loc = purrr::pmap(list(new, wind_tows, sim), ~{
        tow_info <- ..2 |>
          dplyr::filter(sim == ..3, year == 6) |>
          dplyr::select(-c(x, y, cell, depth, AREA_CODE,
                           n, n_aged, n_measured, N))
        dplyr::bind_cols(..1, tow_info)
      })
    ) |>
    dplyr::select(sim, year, strat, new_set_loc) |>
    tidyr::unnest(cols = new_set_loc) |>
    dplyr::relocate(set, .after = dplyr::last_col())

  # Repeat same locations for all post-change years (6–15)
  years_post <- 6:15
  new_locations <- purrr::map_dfr(years_post, function(yr) {
    new_locations_y6 |>
      dplyr::mutate(year = yr)
  })

  # Save output
  saveRDS(
    new_locations,
    here(
      survdat,
      sprintf("%s_%s_%03d_%d_inside_locs_survey.rds",
              species, season, i, nsurveys)
    )
  )
}





#Check
random_locs_001 <- readRDS(here(survdat, "scup_fall_001_25_locs_survey.rds"))
fixed_locs_001 <- readRDS(here(survdat, "scup_fall_001_25_inside_locs_survey.rds"))



# Tows per year
random_locs_001 |> count(year)
fixed_locs_001 |> count(year)

# Tows per year and stratum
random_locs_001 |> count(year, strat)
fixed_locs_001 |> count(year, strat)

# Visualize them
library(ggplot2)
ggplot(random_locs_001, aes(x, y)) +
  geom_point(aes(color = as.factor(year)), alpha = 0.6) +
  facet_wrap(~year) +
  theme_minimal() +
  labs(title = "Reallocated Tow Locations for Pop 001")


ggplot(fixed_locs_001, aes(x, y)) +
  geom_point(aes(color = as.factor(year)), alpha = 0.6) +
  facet_wrap(~year) +
  theme_minimal() +
  labs(title = "Fixed Inside WA Tow Locations for Pop 001")



### LOOP 2: Hybrid survey (random outside wind + reduced fixed inside wind)

set.seed(123)

# define post-change years and fixed sample size
 years_post <- 6:15

 #option 1 for fixed numbers per stratum
 n_fixed_per_stratum <- 1   # number of fixed sites inside wind per stratum

 for (i in this_chunk) {
     message(sprintf("Running hybrid (fixed + random) survey for population %03d of chunk %d...", i, chunk_id))

#Load population and precluded (random outside-WA) survey
pop_i        <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds",species, season, i)))
precl_survey <- readRDS(here(survdat,   sprintf("%s_%s_%03d_%d_precl_survey.rds",species, season, i, nsurveys)))

#Load all inside-wind candidate locations (years 6–15)
inside_all <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_inside_locs_survey.rds",species, season, i, nsurveys)))

#STEP 1: select a smaller fixed subset ONCE (e.g., year 6)
fixed_subset_y6 <- inside_all %>%
   dplyr::filter(AREA_CODE == 1, year == 6) %>%
   dplyr::group_by(sim, strat) %>%
   dplyr::group_modify(~ dplyr::slice_sample(.x,n = min(n_fixed_per_stratum, nrow(.x)),replace = FALSE)) %>%
   dplyr::ungroup()

#STEP 2: repeat the same subset across all post-change years
  fixed_subset_all_years <- fixed_subset_y6 %>%
  dplyr::ungroup() %>%                  # remove hidden grouping
  dplyr::select(-year) %>%              # drop year safely
  tidyr::crossing(year = years_post) %>%
  dplyr::arrange(sim, year, strat) %>%
  dplyr::mutate(set = dplyr::row_number())


#Sanity check
message("  Inside-wind fixed stations: ", nrow(fixed_subset_y6),
   " → repeated across ", length(years_post), " years = ",
 nrow(fixed_subset_all_years), " total rows")

#STEP 3: run SimSurvey for the reduced fixed sites
survey_inside <- sim_survey(
                 pop_i,
                 n_sims          = nsurveys,
                 trawl_dim       = trawl_dim,
                 q               = catch_q,
                 set_den         = set_den,
                 min_sets        = min_sets,
                 age_sampling    = age_sampling,
                 age_space_group = age_space_group,
                 custom_sets     = fixed_subset_all_years,
                 resample_cells  = FALSE)

#STEP 4: combine inside-WA (fixed) and outside-WA (random) surveys
survey_hybrid <- dplyr::bind_rows(precl_survey, survey_inside$setdet)

#STEP 5: save results
saveRDS(fixed_subset_all_years,here(survdat, sprintf("%s_%s_%03d_%d_fixed_locs_reduced_survey.rds",species, season, i, nsurveys)))
saveRDS(survey_hybrid,here(survdat, sprintf("%s_%s_%03d_%d_hybrid_survey.rds",species, season, i, nsurveys)))

message(sprintf("  Saved reduced fixed locations and hybrid survey for population %03d.", i))
                }



 # Example: population 001
 hybrid_001 <- readRDS(here(survdat, "scup_fall_001_25_hybrid_survey.rds"))
 reall_001  <- readRDS(here(survdat, "scup_fall_001_25_reall_survey.rds"))



( compare_sets <- bind_rows(
   hybrid_001  %>% mutate(scenario = "Hybrid"),
   reall_001   %>% mutate(scenario = "Reallocation")
 ) %>%
   count(scenario, year, AREA_CODE) %>%
   pivot_wider(names_from = AREA_CODE,
               values_from = n,
               names_prefix = "area_") %>%
   replace_na(list(area_1 = 0, area_2 = 0)))





 ggplot(hybrid_001, aes(x = x, y = y)) +
       geom_point(aes(color = as.factor(AREA_CODE)), alpha = 0.7) +
       facet_wrap(~year) +
       scale_color_manual(values = c("1" = "#d95f02", "2" = "#1b9e77"),
                          labels = c("Wind area", "Outside wind")) +
       labs(color = "Survey area", title = "Hybrid survey coverage") +
       theme_minimal()


