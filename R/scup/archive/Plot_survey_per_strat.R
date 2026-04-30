##Sandbox for Viz

### DATA SET UP ####
# Directories
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

# Parameters
species   <- "scup"
season    <- "fall"
ages      <- 0:7
years     <- 1:15
nsims     <- 1:100
nsurveys  <- 25
years_post <- 6:15 # Define the years affected by the survey change

#Data
pop <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_abund-dist.rds", species, season, .x))))
survdat_sq <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, .x))))
survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))
survdat_reall <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_reall_survey.rds", species, season, .x))))
dist          <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))

indices <- readRDS(here(surv.prods, str_c(species, season, "all-ihat-25survs-100pops.rds", sep = "_")))


ids     <- sprintf("%03d", nsims)
survdat_hybrid <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_hybrid_survey.rds", species, season, .x))))


# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# find the total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))



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




#indices2 <- bind_rows(ihat_hybrid_all,indices)

tow_summary_hybrid <- map2_dfr(survdat_hybrid, seq_along(survdat_hybrid), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(year, strat) |>
    summarise(
      n_tows = n_distinct(set),
      .groups = "drop"
    ) |>
    mutate(pop = pop_num) |>
    arrange(year, strat)
})

tow_summary_hybrid_wide <- tow_summary_hybrid |>
  pivot_wider(
    names_from = pop,
    values_from = n_tows,
    names_prefix = "pop_"
  ) |>
  arrange(year, strat) |>
  filter (strat == 1050) #filter one strata to check

survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))

tow_summary_precl <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(
      n_tows = n_distinct(set),
      .groups = "drop"
    ) |>
    mutate(pop = pop_num) |>
    arrange(pop, year, strat)
})


tow_summary2_precl <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(year, strat) |>
    summarise(
      n_tows = n_distinct(set),
      .groups = "drop"
    ) |>
    mutate(pop = pop_num) |>
    arrange(year, strat)
})

tow_summary2_precl_wide <- tow_summary2_precl %>%
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

tow_summary_sq2 <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
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

tow_summary_sq_wide2 <- tow_summary_sq2 %>%
  pivot_wider(
    names_from = pop,
    values_from = n_tows,
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



