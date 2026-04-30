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


#BLOCK 1: structured base grid
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

# Add SimSurvey structural attributes (same as in your prior working version)
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

# Now filter INSIDE-WIND cells
grid_inside <- grid_xy |> filter(AREA_CODE == 1)
n_inside <- nrow(grid_inside)


# ---- Match the structure of sim_sets() ----
# sim_sets expands grid × years × sims and includes a unique set ID per row

fixed_inside <- tidyr::expand_grid(
  grid_inside,
  year = years_post,
  sim  = 1:nsurveys) |>
  dplyr::arrange(sim, year, strat, cell) |>
  dplyr::mutate(set = dplyr::row_number()) |>
  dplyr::select(sim, year, strat, x, y, cell, depth, AREA_CODE, division,tow_area, cell_area, strat_cells, strat_area, strat_sets, cell_sets, set)

# ---- Save output ----
saveRDS(
  fixed_inside,
  here("data", "rds", "survdat", "inside_fixed_locs_survey.rds")
)

# Quick structure check
dplyr::glimpse(fixed_inside)


#Check
count <- fixed_inside |>
  group_by(AREA_CODE) |>
  count(year, sim)


# Tows per year
fixed_inside |> count(year)

# Tows per year and stratum
fixed_inside |> count(year, strat, sim)

# Visualize them
ggplot(fixed_inside, aes(x, y)) +
  geom_point(aes(color = as.factor(year)), alpha = 0.6) +
  facet_wrap(~year) +
  theme_minimal() +
  labs(title = "Fixed Inside WA Tow Locations for Pop 001")



grid_counts <- grid_xy %>%
  dplyr::count(AREA_CODE) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      AREA_CODE == 1 ~ "Inside wind area",
      AREA_CODE == 2 ~ "Outside wind area",
      TRUE ~ "Unclassified"
    )
  )




# Summary table of grid coverage by stratum and wind area
strat_summary <- grid_xy %>%
  group_by(strat) %>%
  summarise(
    total_cells = n(),
    wind_cells  = sum(AREA_CODE == 1),  # inside wind areas
    nonwind_cells = sum(AREA_CODE == 2) # outside wind areas
  ) %>%
  ungroup() %>%
  mutate(
    wind_prop = round(wind_cells / total_cells, 3)  # optional: proportion affected
  )


write.csv(strat_summary, here("data", "tables", "stratum_wind_coverage_summary.csv"), row.names = TRUE)


# Count only strata with at least one wind-affected cell
strat_summary_wind <- grid_xy %>%
  group_by(strat) %>%
  summarise(
    total_cells  = n(),
    wind_cells   = sum(AREA_CODE == 1),
    nonwind_cells = sum(AREA_CODE == 2)
  ) %>%
  filter(wind_cells > 0) %>%           # keep only affected strata
  mutate(
    wind_prop = round(wind_cells / total_cells, 3)
  ) %>%
  arrange(desc(wind_prop))             # optional: sort by % affected


# Optional: save table
write.csv(strat_summary_wind,
          here("data", "tables", "stratum_wind_affected_summary.csv"),
          row.names = FALSE)


### BLOCK 2: RUN SUPPLEMENTAL SURVEY (fixed stations inside wind areas) ----
set.seed(123)

# ----- Survey configuration -----
supp_trawl_dim <- c(1.8, 0.010)   # smaller trawl
supp_catch_q   <- sim_logistic(k = 2, x0 = 2.5)


# Define chunking
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))

chunk_id <- 5
this_chunk <- chunks[[chunk_id]]


# Load fixed tow locations (created once from the grid earlier)
fixed_inside <- readRDS(here("data", "rds", "survdat", "inside_fixed_locs_survey.rds"))

# ----- Loop over populations in this chunk -----
for (i in this_chunk) {
  message(sprintf("Running supplemental (fixed) survey for population %03d of chunk %d...", i, chunk_id))

  # Load population abundance-distribution object
  pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds",
                                          species, season, i)))

  # Run survey for fixed stations only (no random sampling, no reallocation)
  survey_fixed <- sim_survey(
    sim              = pop_i,
    n_sims           = nsurveys,                 # 25 supplemental survey replicates
    trawl_dim        = supp_trawl_dim,
    q                = supp_catch_q,
    resample_cells   = FALSE,              # fixed locations only
    custom_sets      = fixed_inside,       # use your inside-wind grid
    age_sampling     = "random",
    age_space_group  = "set"
  )

  # Save result
  saveRDS(
    survey_fixed,
    here(survdat, sprintf("%s_%s_%03d_%d_supplemental_survey.rds",
                          species, season, i, nsurveys))
  )

  message(sprintf("  Saved supplemental survey for population %03d.", i))
}




 # Example: population 001
 suppl_001 <- readRDS(here(survdat, "scup_fall_001_25_supplemental_survey.rds"))



( sets <- bind_rows(
  suppl_001$setdet  |> mutate(scenario = "Supplemental")) |>
   count(scenario, year, AREA_CODE, sim) %>%
   pivot_wider(names_from = AREA_CODE,
               values_from = n,
               names_prefix = "area_") %>%
   replace_na(list(area_1 = 0, area_2 = 0)))



 # choose a few populations to inspect
 pops_to_check <- c(1, 5, 12)

 survdat <- here("data", "rds", "survdat")

 supplemental_list <- map(pops_to_check, function(i) {
   readRDS(here(survdat, sprintf("scup_fall_%03d_25_supplemental_survey.rds", i)))
 })
 names(supplemental_list) <- sprintf("Pop_%03d", pops_to_check)


 bind_rows(
   lapply(names(supplemental_list), function(nm) {
     supplemental_list[[nm]]$setdet %>%
       mutate(pop = nm)
   })
 ) %>%
   ggplot(aes(x = x, y = y, color = pop)) +
   geom_point(alpha = 0.6, size = 1) +
   facet_grid(sim~pop) +
   theme_bw() +
   labs(title = "Supplemental Surveys Across Subset Populations",
        subtitle = "") +
   theme(
     legend.position = "none",
     strip.text = element_text(size = 12, face = "bold"),
     axis.title = element_blank()
   )



 supp_abund <- bind_rows(
   lapply(names(supplemental_list), function(nm) {
     supplemental_list[[nm]]$setdet %>%
       group_by(year, sim) %>%
       summarise(total_biomass = sum(N, na.rm = TRUE)) %>%
       mutate(pop = nm)
   })
 )


 ggplot(supp_abund, aes(x = year, y = total_biomass, color = pop, group = pop)) +
   geom_line(linewidth = .8) +
   geom_point(size = 2) +
   facet_wrap(~pop, scales = "free_y") +
   theme_bw() +
   theme(legend.position = "none") +
   labs(
     title = "Supplemental Survey Abundance (Populations 001, 005, 012)",
     x = "Year",
     y = "Total Biomass (kg)"
   )




 ### BLOCK 3: Run HYBRID survey (preclusion + supplemental)
 chunk_id <- 5
 this_chunk <- chunks[[chunk_id]]


 for (i in this_chunk) {
   message(sprintf("Running HYBRID (preclusion + supplemental) survey for population %03d of chunk %d...", i, chunk_id))

   #Load objects ---
   precl_survey <- readRDS(here(survdat,sprintf("%s_%s_%03d_%d_precl_survey.rds",species, season, i, nsurveys)))
   suppl_survey <- readRDS(here(survdat,sprintf("%s_%s_%03d_%d_supplemental_survey.rds",species, season, i, nsurveys)))

   # Combine both survey outputs, Keep only the setdet table from each (SimSurvey outputs a list)
   survey_hybrid <- dplyr::bind_rows(precl_survey, suppl_survey$setdet)

   # --- Save the hybrid survey ---
   saveRDS(survey_hybrid,here(survdat,sprintf("%s_%s_%03d_%d_hybrid_survey.rds",species, season, i, nsurveys)))

   message(sprintf("Saved HYBRID survey for population %03d (%d rows total)",i, nrow(survey_hybrid)))
 }


 ids <- c(1, 5, 12)
  hybrid_list <- setNames(
  lapply(ids, function(i) {
  readRDS(here("data", "rds", "survdat", sprintf("scup_fall_%03d_25_hybrid_survey.rds", i)))
         }),
        sprintf("Pop_%03d", ids)
    )

 hybrid_all <- data.table::rbindlist(
      lapply(names(hybrid_list), function(nm) {
          df <- as.data.frame(hybrid_list[[nm]])
          df$pop <- nm
          return(df)
      }),
      use.names = TRUE, fill = TRUE)




 ggplot(hybrid_all, aes(x = x, y = y)) +
   geom_point(aes(color = as.factor(AREA_CODE)), alpha = 0.7, size = 0.8) +
   facet_wrap(~pop) +
   scale_color_manual(values = c("1" = "#d95f02", "2" = "#1b9e77"),
                      labels = c("Wind area", "Outside wind")) +
   coord_equal() +
   labs(
     title = "Hybrid Survey Subset",
     x = "Longitude", y = "Latitude",
     color = "Survey area"
   ) +
   theme_minimal(base_size = 13) +
   theme(legend.position = "bottom")



 ### DATA SET UP ####
 # Directories
 sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
 dist.dat <- here("data", "rds", "dists")
 survdat <- here("data", "rds", "survdat")
 surv.prods <- here("data", "rds", "surv-prods")
 plots <- here("outputs", "plots")

 ids     <- sprintf("%03d", nsims)
 survdat_hybrid <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_hybrid_survey.rds", species, season, .x))))


 # area weights for each strata
 strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
   rename(strat = STRATUM)

 # find the total survey area
 survey_area <- as.integer(sum(strata_wts$Area_SqNm))


 #function to set selectivity
 source(here("R/selectivity_fns.R"))
 q = force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1)
 (selectivity_values <- q(ages))


 ### Hybrid Survey ####
 #FOR ALL THE REALIZATIONS OF THE POPULATION
 ihat_hybrid <- map2_dfr(survdat_hybrid, seq_along(survdat_hybrid), function(surv, pop_num) {
   surv |> #it is already a data.table
     as_tibble() |>
     filter(strat %in% unlist(strat)) |>
     group_by(sim, year, strat) |>
     summarise(towct = length(unique(set)),
               mu = sum(n)/towct,
               var = ifelse(towct == 1, 0,
                            sum((n - mu)^2)/(towct - 1)),
               .groups = "drop") |>
     left_join(strata_wts, by = "strat") |>
     mutate(wt_mu = Area_SqNm * mu,
            wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
     group_by(sim, year) |>
     summarise(stratmu = (sum(wt_mu)) / survey_area,
               stratvar = sum(wt_var),
               cv = sqrt(stratvar)/stratmu,
               .groups = "drop") |>
     mutate(scenario = "Hybrid",
            pop = pop_num)
 })

 # Group by pop (and sim) to compute rel_ihat within each population
 ihat_hybrid_all <- ihat_hybrid %>%
   group_by(pop, sim) %>%
   mutate(rel_ihat = stratmu / mean(stratmu),
          log_se = sqrt(log(1 + cv^2)),
          log_lower = log(stratmu) - 1.96 * log_se,
          log_upper = log(stratmu) + 1.96 * log_se,
          ci_lower = exp(log_lower),
          ci_upper = exp(log_upper)) %>%
   ungroup()







 # Assumes: hybrid_all is your merged data.table of hybrid surveys
 #          strata_wts has columns: strat, Area_SqNm, RelWt, survey_area



 indices2 <- bind_rows(ihat_hybrid,indices)



 tow_summary_hybrid <- map2_dfr(survdat_hybrid, seq_along(survdat_hybrid), function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat %in% unlist(strat)) %>%
     group_by(year, strat) %>%
     summarise(
       n_tows = n_distinct(set),
       .groups = "drop"
     ) %>%
     mutate(pop = pop_num) %>%
     arrange(year, strat)
 })

 tow_summary_hybrid_wide <- tow_summary_hybrid %>%
   pivot_wider(
     names_from = pop,
     values_from = n_tows,
     names_prefix = "pop_"
   ) %>%
   arrange(year, strat) |>
   filter (strat == 1050)


 tow_summary_precl <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat %in% unlist(strat)) %>%
     group_by(sim, year, strat) %>%
     summarise(
       n_tows = n_distinct(set),
       .groups = "drop"
     ) %>%
     mutate(pop = pop_num) %>%
     arrange(pop, year, strat)
 })


 tow_summary_precl <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat %in% unlist(strat)) %>%
     group_by(year, strat) %>%
     summarise(
       n_tows = n_distinct(set),
       .groups = "drop"
     ) %>%
     mutate(pop = pop_num) %>%
     arrange(year, strat)
 })

 tow_summary_precl_wide <- tow_summary_precl %>%
   pivot_wider(
     names_from = pop,
     values_from = n_tows,
     names_prefix = "pop_"
   ) %>%
   arrange(year, strat) |>
   filter (strat == 1050)


 tow_summary_reall <- map2_dfr(survdat_reall, seq_along(survdat_reall), function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat %in% unlist(strat)) %>%
     group_by(sim, year, strat) %>%
     summarise(
       n_tows = n_distinct(set),
       .groups = "drop"
     ) %>%
     mutate(pop = pop_num) %>%
     arrange(pop, year, strat)
 })


 tow_summary_sq <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
   surv$setdet %>%
     as_tibble() %>%
     filter(strat %in% unlist(strat)) %>%
     group_by(sim, year, strat) %>%
     summarise(
       n_tows = n_distinct(set),
       .groups = "drop"
     ) %>%
     mutate(pop = pop_num) %>%
     arrange(pop, year, strat)
 })

 tow_summary_sq <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
   surv$setdet %>%
     as_tibble() %>%
     filter(strat %in% unlist(strat)) %>%
     group_by(year, strat) %>%
     summarise(
       n_tows = n_distinct(set),
       .groups = "drop"
     ) %>%
     mutate(pop = pop_num) %>%
     arrange(year, strat)
 })

 tow_summary_sq_wide2 <- tow_summary_sq_2 %>%
   pivot_wider(
     names_from = pop,
     values_from = n_sets,
     names_prefix = "pop_"
   ) %>%
   arrange(year, strat) |>
   filter (strat == 1050)







 # read the strata shapefile in to plot the polygons based on proportion and overlay the tow points
 strata <- readRDS(here(sseep.analysis,"data", "rds", "active_strata.rds"))
 wind_areas <- readRDS(here(sseep.analysis,"data", "rds", "all_wind_areas_Jun2022.rds"))
 strat_1050 <- strata %>%
   filter(STRATUM == 1050) %>%
   st_transform(crs = 4326)   # convert to same CRS as your other layers

 # Pick the stratum and 5 populations you want to visualize
 target_strat <- 1050
 target_pops  <- 1:3  # you can change to specific pop IDs

 # Combine the $setdet data for the selected populations
 tow_subset <- map2_dfr(survdat_sq[target_pops], target_pops, function(surv, pop_num) {
   surv$setdet %>%
     as_tibble() %>%
     filter(strat == target_strat) %>%    # select the stratum of interest
     mutate(pop = pop_num)
 })


 tow_subset <- tow_subset |> filter(sim %in% c(3,6,9),  year %in% c(1,5,6,10,15)) |>
   mutate(x_m = x * 1000, y_m = y * 1000)

 # Correct UTM zone for the U.S. Northeast: 18N
 tow_subset_sf <- st_as_sf(tow_subset, coords = c("x_m", "y_m"), crs = 32618) |>
   st_transform(crs = 4326)

#st_bbox(tow_subset_sf)
 ggplot() +
   geom_sf(data = wind_areas, fill = "white", color = "grey20", alpha = 0.5) +
   geom_sf(data = tow_subset_sf, aes(color = factor(pop)), size = 1.5, alpha = 0.8) +
   geom_sf(data = strat_1050, fill = NA, color = "firebrick2", linewidth = 1) +
   facet_grid(sim~year) +
   coord_sf(xlim = c(-73, -70), ylim = c(40, 42), expand = FALSE) +
   labs(
     title = paste("Tow locations in stratum", target_strat),
     subtitle = "Status Quo Survey",
     x = "Longitude", y = "Latitude",
     color = "Population"
   ) +
   theme_bw()






 # Combine the $setdet data for the selected populations
 tow_subset_precl <- map2_dfr(survdat_precl[target_pops], target_pops, function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat == target_strat) %>%    # select the stratum of interest
     mutate(pop = pop_num)
 })


 tow_subset_precl <- tow_subset_precl |> filter(sim %in% c(3,6,9),  year %in% c(1,5,6,10,15)) |>
   mutate(x_m = x * 1000, y_m = y * 1000)

 # Correct UTM zone for the U.S. Northeast: 18N
 tow_subset_precl_sf <- st_as_sf(tow_subset_precl, coords = c("x_m", "y_m"), crs = 32618) |>
   st_transform(crs = 4326)

 #st_bbox(tow_subset_sf)
 ggplot() +
   geom_sf(data = wind_areas, fill = "white", color = "grey20", alpha = 0.5) +
   geom_sf(data = tow_subset_precl_sf, aes(color = factor(pop)), size = 1.5, alpha = 0.8) +
   geom_sf(data = strat_1050, fill = NA, color = "firebrick2", linewidth = 1) +
   facet_grid(sim~year) +
   coord_sf(xlim = c(-73, -70), ylim = c(40, 42), expand = FALSE) +
   labs(
     title = paste("Tow locations in stratum", target_strat),
     subtitle = "Preclusion Survey",
     x = "Longitude", y = "Latitude",
     color = "Population"
   ) +
   theme_bw()






 # Combine the $setdet data for the selected populations
 tow_subset_reall <- map2_dfr(survdat_reall[target_pops], target_pops, function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat == target_strat) %>%    # select the stratum of interest
     mutate(pop = pop_num)
 })


 tow_subset_reall <- tow_subset_reall |> filter(sim %in% c(3,6,9),  year %in% c(1,5,6,10,15)) |>
   mutate(x_m = x * 1000, y_m = y * 1000)

 # Correct UTM zone for the U.S. Northeast: 18N
 tow_subset_reall_sf <- st_as_sf(tow_subset_reall, coords = c("x_m", "y_m"), crs = 32618) |>
   st_transform(crs = 4326)

 #st_bbox(tow_subset_sf)
 ggplot() +
   geom_sf(data = wind_areas, fill = "white", color = "grey20", alpha = 0.5) +
   geom_sf(data = tow_subset_reall_sf, aes(color = factor(pop)), size = 1.5, alpha = 0.8) +
   geom_sf(data = strat_1050, fill = NA, color = "firebrick2", linewidth = 1) +
   facet_grid(sim~year) +
   coord_sf(xlim = c(-73, -70), ylim = c(40, 42), expand = FALSE) +
   labs(
     title = paste("Tow locations in stratum", target_strat),
     subtitle = "Reallocation Survey",
     x = "Longitude", y = "Latitude",
     color = "Population"
   ) +
   theme_bw()





  # Combine the $setdet data for the selected populations
 tow_subset_hybrid <- map2_dfr(survdat_hybrid[target_pops], target_pops, function(surv, pop_num) {
   surv %>%
     as_tibble() %>%
     filter(strat == target_strat) %>%    # select the stratum of interest
     mutate(pop = pop_num)
 })


 tow_subset_hybrid <- tow_subset_hybrid |> filter(sim %in% c(3,6,9),  year %in% c(1,5,6,10,15)) |>
   mutate(x_m = x * 1000, y_m = y * 1000)

 # Correct UTM zone for the U.S. Northeast: 18N
 tow_subset_hybrid_sf <- st_as_sf(tow_subset_hybrid, coords = c("x_m", "y_m"), crs = 32618) |>
   st_transform(crs = 4326)

 #st_bbox(tow_subset_sf)
 ggplot() +
   geom_sf(data = wind_areas, fill = "white", color = "grey20", alpha = 0.5) +
   geom_sf(data = tow_subset_hybrid_sf, aes(color = factor(pop)), size = 1.5, alpha = 0.8) +
   geom_sf(data = strat_1050, fill = NA, color = "firebrick2", linewidth = 1) +
   facet_grid(sim~year) +
   coord_sf(xlim = c(-73, -70), ylim = c(40, 42), expand = FALSE) +
   labs(
     title = paste("Tow locations in stratum", target_strat),
     subtitle = "Hybrid Survey",
     x = "Longitude", y = "Latitude",
     color = "Population"
   ) +
   theme_bw()



