### created: 10/14/2025
### updated: 10/30/2025

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
library(sf)
library(dplyr)
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
years_post <- 6:15 # Define the years affected by the survey change


###BLOCK 1: structured base grid
grid_xy <- readRDS(here("data", "rds", "survey_grid_all_122024.rds")) |>
  rename(
    strat = STRATUM,
    x = X,
    y = Y,
    depth = mean_2
  ) |>
  mutate(division = 1) |>
  dplyr::select(x, y, cell, depth, strat, AREA_CODE, division) |>
  data.table::as.data.table() |>
  drop_na()

# a) Add SimSurvey structural attributes

set_den <- 0.001

grid_xy <- grid_xy |>
  group_by(strat) |>
  mutate(
    tow_area    = 0.0378,
    cell_area   = 16.2,
    strat_cells = n(),
    strat_area  = strat_cells * cell_area,
    strat_sets  = round(strat_area * set_den),
    strat_sets  = ifelse(strat_sets < 3, 3, strat_sets),
    cell_sets   = 1
  ) |>
  ungroup()

# b) filter INSIDE-WIND cells
grid_inside <- grid_xy |> filter(AREA_CODE == 1)
n_inside <- nrow(grid_inside)


# c) Match the structure of sim_sets()
# sim_sets expands grid × years × sims and includes a unique set ID per row

fixed_inside <- tidyr::expand_grid(
  grid_inside,
  year = years_post,
  sim  = 1:nsurveys) |>
  dplyr::arrange(sim, year, strat, cell) |>
  dplyr::mutate(set = dplyr::row_number()) |>
  dplyr::select(sim, year, strat, x, y, cell, depth, AREA_CODE, division,tow_area, cell_area, strat_cells, strat_area, strat_sets, cell_sets, set)


# d) Save output
saveRDS(fixed_inside, here("data", "rds", "survdat", "inside_fixed_locs_survey.rds"))


#Optional just for check structure, summary and viz
# #Check
# # structure check
# dplyr::glimpse(fixed_inside)
#
# count <- fixed_inside |>
#   group_by(AREA_CODE) |>
#   count(year, sim)
#
# # Tows per year
# fixed_inside |> count(year)
#
# # Tows per year and stratum
# fixed_inside |> count(year, strat, sim)
#
# # Visualize them
# ggplot(fixed_inside, aes(x, y)) +
#   geom_point(aes(color = as.factor(year)), alpha = 0.6) +
#   facet_wrap(~year) +
#   theme_minimal() +
#   labs(title = "Fixed Inside WA Tow Locations for Pop 001")
#
#
# grid_counts <- grid_xy |>
#   dplyr::count(AREA_CODE) |>
#   dplyr::mutate(
#     category = dplyr::case_when(
#       AREA_CODE == 1 ~ "Inside wind area",
#       AREA_CODE == 2 ~ "Outside wind area",
#       TRUE ~ "Unclassified"))
#
# # Summary table of grid coverage by stratum and wind area
# strat_summary <- grid_xy |>
#   group_by(strat) |>
#   summarise(
#     total_cells = n(),
#     wind_cells  = sum(AREA_CODE == 1),  # inside wind areas
#     nonwind_cells = sum(AREA_CODE == 2) # outside wind areas
#   ) |>
#   ungroup() |>
#   mutate(wind_prop = round(wind_cells / total_cells, 3))  # proportion affected
#
# write.csv(strat_summary, here("data", "tables", "stratum_wind_coverage_summary.csv"), row.names = TRUE)
#
#
# # Count only strata with at least one wind-affected cell
# strat_summary_wind <- grid_xy |>
#   group_by(strat) |>
#   summarise(
#     total_cells  = n(),
#     wind_cells   = sum(AREA_CODE == 1),
#     nonwind_cells = sum(AREA_CODE == 2)
#   ) |>
#   filter(wind_cells > 0)|>         # keep only affected strata
#   mutate(wind_prop = round(wind_cells / total_cells, 3)) |>
#   arrange(desc(wind_prop))             # optional: sort by % affected
#
#
# # save table
# write.csv(strat_summary_wind, here("data", "tables", "stratum_wind_affected_summary.csv"), row.names = TRUE)



### BLOCK 2: RUN SUPPLEMENTAL SURVEY (fixed stations inside wind areas)

# a) Survey configuration
supp_trawl_dim <- c(1.8, 0.010)   # smaller trawl
supp_catch_q   <- sim_logistic(k = 2, x0 = 2.5) #logistic, different as used for scup in SQ survey


# b) Define chunking
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))

chunk_id <- 1 #change from 1 to 5
this_chunk <- chunks[[chunk_id]]


# c) Load fixed tow locations (created here in this script)
fixed_inside <- readRDS(here("data", "rds", "survdat", "inside_fixed_locs_survey.rds"))

# d) Loop over populations in this chunk
set.seed(123)
for (i in this_chunk) {
  message(sprintf("Running supplemental (fixed) survey for population %03d of chunk %d...", i, chunk_id))

  # Load population abundance-distribution object
  pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

  # Run survey for fixed stations only (no random sampling, no reallocation)
  survey_fixed <- sim_survey(
    sim              = pop_i,
    n_sims           = nsurveys,           # 25 identical supplemental survey replicates
    trawl_dim        = supp_trawl_dim,
    q                = supp_catch_q,
    resample_cells   = FALSE,              # fixed locations only
    custom_sets      = fixed_inside,       # use your inside-wind grid
    age_sampling     = "random",
    age_space_group  = "set"
  )

  # Save result
  saveRDS(survey_fixed, here(survdat, sprintf("%s_%s_%03d_%d_supplemental_survey.rds", species, season, i, nsurveys)))
  message(sprintf("  Saved supplemental survey for population %03d.", i))
}



# Optional check + viz
# Example: population 001
# suppl_001 <- readRDS(here(survdat, "scup_fall_001_25_supplemental_survey.rds"))
#
#
# (sets <- bind_rows(
#   suppl_001$setdet  |> mutate(scenario = "Supplemental")) |>
#    count(scenario, year, AREA_CODE, sim) |>
#    pivot_wider(names_from = AREA_CODE,
#                values_from = n,
#                names_prefix = "area_") |>
#    replace_na(list(area_1 = 0, area_2 = 0)))
#
#
#
#  # choose a few populations to inspect
#  pops_to_check <- c(1, 5, 12)
#
#  supplemental_list <- map(pops_to_check, function(i) {
#    readRDS(here(survdat, sprintf("scup_fall_%03d_25_supplemental_survey.rds", i)))
#  })
#
#  names(supplemental_list) <- sprintf("Pop_%03d", pops_to_check)
#
#
#  bind_rows(
#    lapply(names(supplemental_list), function(nm) {
#      supplemental_list[[nm]]$setdet |>
#        mutate(pop = nm)
#    })
#  ) |>
#    ggplot(aes(x = x, y = y, color = pop)) +
#    geom_point(alpha = 0.6, size = 1) +
#    facet_grid(sim~pop) +
#    theme_bw() +
#    labs(title = "Supplemental Surveys Across Subset Populations",
#         subtitle = "") +
#    theme(
#      legend.position = "none",
#      strip.text = element_text(size = 12, face = "bold"),
#      axis.title = element_blank()
#    )
#
#
#  supp_abund <- bind_rows(
#    lapply(names(supplemental_list), function(nm) {
#      supplemental_list[[nm]]$setdet |>
#        group_by(year, sim) |>
#        summarise(total_biomass = sum(N, na.rm = TRUE)) |>
#        mutate(pop = nm)
#    })
#  )
#
#
#  ggplot(supp_abund, aes(x = year, y = total_biomass, color = pop, group = pop)) +
#    geom_line(linewidth = .8) +
#    geom_point(size = 2) +
#    facet_wrap(~pop, scales = "free_y") +
#    theme_bw() +
#    theme(legend.position = "none") +
#    labs(
#      title = "Supplemental Survey Abundance (Populations 001, 005, 012)",
#      x = "Year",
#      y = "Total Biomass (kg)"
#    )




### BLOCK 3: Arrange HYBRID survey (preclusion + supplemental)
# a) Define chunking
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))
chunk_id <- 1 #change from 1 to 5
this_chunk <- chunks[[chunk_id]]


# b) Combine preclusion (years 1-5 all and outside WA) + fixed survey (years 6-15 inside WA)
 for (i in this_chunk) {
   message(sprintf("Running HYBRID (preclusion + supplemental) survey for population %03d of chunk %d...", i, chunk_id))

   #Load objects ---
   precl_survey <- readRDS(here(survdat,sprintf("%s_%s_%03d_%d_precl_survey.rds",species, season, i, nsurveys)))
   suppl_survey <- readRDS(here(survdat,sprintf("%s_%s_%03d_%d_supplemental_survey.rds",species, season, i, nsurveys)))

   # Combine both survey outputs, Keep only the setdet table from each
   survey_hybrid <- dplyr::bind_rows(precl_survey, suppl_survey$setdet)

   # Save
   saveRDS(survey_hybrid,here(survdat,sprintf("%s_%s_%03d_%d_hybrid_survey.rds",species, season, i, nsurveys)))

   message(sprintf("Saved HYBRID survey for population %03d (%d rows total)",i, nrow(survey_hybrid)))
 }



# Optional check + viz
# ##make sure that preclusion survey is being added properly to the "hybrid" survey
# idss <- c(25,75)
#
# # load the preclusion surveys
# precl_list <- setNames(
#  lapply(idss, \(i)
#  readRDS(here("data", "rds", "survdat", sprintf("scup_fall_%03d_25_precl_survey.rds", i)))),
#  sprintf("Pop_%03d", idss))
#
# # load the hybrid surveys
# hybrid_list <- setNames(
#  lapply(idss, \(i)
#  readRDS(here("data", "rds", "survdat", sprintf("scup_fall_%03d_25_hybrid_survey.rds", i)))),
#  sprintf("Pop_%03d", idss))
#
#  summarise_survey <- function(df) {
#    df |>
#      filter(year == 10) |>
#      group_by(year) |>
#      summarise(
#        n_tows = n_distinct(set),
#        total_N = sum(N, na.rm = TRUE),
#        .groups = "drop"
#      )}
#
#  precl_summary  <- map_df(precl_list, summarise_survey, .id = "Population") |>
#    mutate(Survey = "Preclusion")
#
#  hybrid_summary <- map_df(hybrid_list, summarise_survey, .id = "Population") |>
#    mutate(Survey = "Hybrid")
#
#
# ids <- c(1, 5, 12) #population subset
# hybrid_list <- setNames(
# lapply(ids, function(i) {
# readRDS(here("data", "rds", "survdat", sprintf("scup_fall_%03d_25_hybrid_survey.rds", i)))
#   }),
#   sprintf("Pop_%03d", ids))
#
# hybrid_all <- data.table::rbindlist(
#     lapply(names(hybrid_list), function(nm) {
#         df <- as.data.frame(hybrid_list[[nm]])
#         df$pop <- nm
#         return(df)
#       }),
#       use.names = TRUE, fill = TRUE)
#
#
# ggplot(hybrid_all, aes(x = x, y = y)) +
#    geom_point(aes(color = as.factor(AREA_CODE)), alpha = 0.7, size = 0.8) +
#    facet_wrap(~pop) +
#    scale_color_manual(values = c("1" = "#d95f02", "2" = "#1b9e77"),
#                       labels = c("Wind area", "Outside wind")) +
#    coord_equal() +
#    labs(
#      title = "Hybrid Survey Subset",
#      x = "Longitude", y = "Latitude",
#      color = "Survey area"
#    ) +
#    theme_minimal(base_size = 13) +
#    theme(legend.position = "bottom")



