### created: 05/16/2025
### updated:

# 04c - SIMULATE SUPPLEMENTAL SAMPLING ####
# random stratified outside
# fixed stations inside wind areas


## Objective ####


# Outputs:

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
#chunk_size <- 20
#chunks <- split(nsims, ceiling(nsims / chunk_size))


# Trawl survey configuration
trawl_dim       <- c(2.7, 0.014)
catch_q         <- force_sim_logistic(k = -0.66, x0 = -1.14, force_age = TRUE, age = 0, force_sel = 1)
set_den         <- 0.001
min_sets        <- 3
age_sampling    <- "stratified"
age_space_group <- "set"
resample_cells  <- TRUE


# Set which chunk to run (1 to 5)
#chunk_id <- Change this (1 to 5)
#this_chunk <-  chunks[[chunk_id]]

# Parameters
species   <- "scup"
season    <- "fall"
ages      <- 0:7
years     <- 1:15
nsims     <- 1:100
nsurveys  <- 25
this_chunk <- list(1)

#LOOP 1 Generate 25 different sets of fixed tow locations inside wind areas

for (i in this_chunk) {
  message(sprintf("Generating fixed wind-area locations for population %03d...", i))

  # Load SQ survey
  sq_survey <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_sq_survey.rds", species, season, i, nsurveys)))
  grid      <- sq_survey$grid_xy
  setdet    <- sq_survey$setdet

  # Identify wind tows after year 5 and extract strata
  wind_tows   <- setdet |> filter(AREA_CODE == 1, year >= 6)
  wind_strata <- unique(wind_tows$strat)

  # Filter grid INSIDE wind areas for affected strata
  grid_wind <- grid |> filter(AREA_CODE == 1, strat %in% wind_strata)
  in_wa_grid <- grid_wind |> group_by(strat) |> nest()
  in_strat <- unique(in_wa_grid$strat)

  # Count wind tows per sim-year-stratum
  wind_summ <- wind_tows |>
    group_by(sim, year, strat) |>
    nest() |>
    mutate(count = map(data, ~length(.$set))) |>
    rename(wind_tows = data)

  # Join with available wind grid cells and keep strata with valid area
  join_data <- left_join(wind_summ, in_wa_grid, by = "strat") |>
    filter(strat %in% in_strat) |>
    filter(!map_lgl(data, is.null))

  # Generate fixed locations INSIDE wind areas
  new_locations <- join_data |>
    mutate(
      new = map2(data, count, ~slice_sample(.x, n = .y, replace = TRUE)),
      new_set_loc = pmap(list(new, wind_tows, sim), ~{
        tow_info <- ..2 |> filter(sim == ..3) |>
          select(-c(x, y, cell, depth, AREA_CODE, n, n_aged, n_measured, N))
        bind_cols(..1, tow_info)
      })
    ) |>
    select(sim, year, strat, new_set_loc) |>
    unnest(cols = new_set_loc) |>
    relocate(set, .after = last_col())

  # Save output
  saveRDS(new_locations, here(survdat, sprintf("%s_%s_%03d_%d_locs_supp_wind_fixed.rds", species, season, i, nsurveys)))
}




#Check
supp_locs <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_locs_supp_wind_fixed.rds", species, season, 1, nsurveys)))
#glimpse(supp_locs)

supp_locs |>
  group_by(year) |>
  summarise(n_unique_cells = n_distinct(cell))



# Convert to sf in meters
supp_locs_sf <- supp_locs |>
  mutate(x = x * 1000, y = y * 1000) |>
  st_as_sf(coords = c("x", "y"), crs = 32618)


sseep.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
strata <- readRDS(here(sseep.dir,"data", "rds", "active_strata.rds"))
wind_areas <- readRDS(here(sseep.dir,"data", "rds", "all_wind_areas_Jun2022.rds"))

strata_proj <- st_transform(strata, crs = 32618)
wind_proj   <- st_transform(wind_areas, crs = 32618)
library(sf)

fivesims <- supp_locs_sf |> filter(sim %in% 1:5)


ggplot() +
  geom_sf(data = wind_proj, fill = "lightblue", alpha = 0.4, color = NA) +
  geom_sf(data = strata_proj, fill = NA, color = "gray50", linewidth = 0.3) +
  geom_sf(data = fivesims, color = "orange", size = 1) +
  facet_wrap(~sim) +
  coord_sf(crs = 32618, datum = NA, expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Supplemental Tow Locations by Year",
    subtitle = "Fixed supplemental stations inside wind areas (UTM Zone 18N)",
    x = "Easting (m)", y = "Northing (m)"
  )



# Tows per year
supp_locs |> count(year, sim)

# Tows per year and stratum
supp_locs |> count(year, strat, sim)

# Count number of tows per strat, year, sim
tow_counts <- supp_locs_sf |>
  st_drop_geometry() |>  # 👈 DROP geometry first!
  count(strat, year, sim, name = "n_tows")

# Reshape to wide format: one row per strat+year, one column per sim
tow_table <- tow_counts |>
  pivot_wider(names_from = sim, values_from = n_tows, values_fill = 0)




#LOOP 2: un Survey at Fixed Supplemental Locations
# Parameters
this_chunk <- list(1)
i=1

for (i in this_chunk) {
  message(sprintf("Running supplemental fixed-location survey for population %03d...", i))

  # Load population and fixed tow locations
  pop_i     <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))
  supp_locs <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_locs_supp_wind_fixed.rds", species, season, i, nsurveys)))

  # Run survey with those fixed tow locations
  supp_survey <- sim_survey(
    sim               = pop_i,
    n_sims            = nsurveys,
    trawl_dim         = trawl_dim,
    q                 = catch_q,
    set_den           = set_den,
    min_sets          = min_sets,
    age_sampling      = age_sampling,
    age_space_group   = age_space_group,
    custom_sets       = supp_locs,
    resample_cells    = resample_cells
  )

  # Save just the fixed-location supplemental survey
  saveRDS(supp_survey$setdet, here(survdat, sprintf("%s_%s_%03d_%d_supp_wind_only_survey.rds", species, season, i, nsurveys)))

  # Load precluded survey (from SQ - wind)
  precl_survey <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_precl_survey.rds", species, season, i, nsurveys)))

  # Merge both surveys (same structure $setdet)
  survey_supp_combined <- bind_rows(precl_survey, supp_survey$setdet)

  # Save merged result
  saveRDS(survey_supp_combined, here(survdat, sprintf("%s_%s_%03d_%d_supp_wind_fixed_survey.rds", species, season, i, nsurveys)))
}


survdat_fixed_1 <- readRDS(here(survdat, "scup_fall_001_25_supp_wind_only_survey.rds"))

survdat_comb_1<- readRDS(here(survdat, "scup_fall_001_25_supp_wind_fixed_survey.rds"))

survdat_reall_1 <- readRDS(here(survdat, "scup_fall_001_25_reall_survey.rds"))





survdat_reall_old[[1]]




#Count
s_comb <- survdat_comb_1

# Count sets by sim, year
tow_counts3 <- s_comb |>
  count(sim, year, name = "n_sets") #can add strat to check

print(tow_counts3)


#Average number of tows per year across all simulations
avg_tows_per_year3 <- s_comb |>
  count(sim, year) |>
  group_by(year) |>
  summarise(
    mean_tows = mean(n),
    sims = n(),
    .groups = "drop"
  )

print(avg_tows_per_year3)












# Get wind areas intersecting stratum 3380
strat_3380 <- strata_proj |> filter(STRATUM == 3380)
wind_areas_3380 <- wind_proj[st_intersects(wind_proj, strat_3380, sparse = FALSE)[,1], ]
# Get full geometry of all strata that intersect those wind areas
strata_touching_same_wind <- strata_proj[st_intersects(strata_proj, wind_areas_3380, sparse = FALSE) |> apply(1, any), ]
# Plot background for context (optional coastline or bounding box)
plot(st_geometry(strata_touching_same_wind), col = "lightblue", border = "blue", main = "Strata 3380 and 1650 Intersecting Wind Area")
# Add wind areas
plot(st_geometry(wind_areas_3380), col = adjustcolor("red", alpha.f = 0.5), border = "darkred", add = TRUE)
# Add stratum 3380 with thicker orange outline
plot(st_geometry(strat_3380), border = "orange", col = NA, lwd = 2, add = TRUE)
strata_touching_same_wind




















