### created: 01/18/2024
### updated: 02/07/2024

# 05 - CALCULATE RELATIVE ABUNDANCES ####


## Objective ####
# For a given species, distribution, and survey, calculate the relative true abundance and the relative abundance index.

# Outputs: Relative true abundance and abundance indices for each year and simulation of the projection standardized to the average abundance or abundance index over time

### PACKAGES ####
library(tidyverse)
library(here)
library(sdmTMB)
library(SimSurvey)
source(here("R", "sim_stratmean_fn.R"))



### DATA SET UP ####
# Directories
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

# Parameters
species <- "summerflounder"
season  <- "fall"
ages      <- 0:7
years     <- 1:15
nsims   <- 1:100
ids     <- sprintf("%03d", nsims)

#Data
pop <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_abund-dist.rds", species, season, .x))))
survdat_sq <- map(ids, function(id) {x <- readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, id)))
  out <- x$setdet #load only setdet data
  rm(x); gc()
  out
})
survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))
survdat_reall <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_reall_survey.rds", species, season, .x))))
dist          <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))




# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# find the total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))


#function to set selectivity
source(here("R/selectivity_fns.R"))

q <- sim_logistic(k = 2, x0 = 2.5)
ages <- 0:7
selectivity_values <- q(ages)
names(selectivity_values) <- ages
selectivity_values


## True Abundance ####
# calculate the relative true abundance from the simulated population and distribution
trueN <- map(pop, ~as_tibble(.$N) |>
               mutate(age = ages) |>
               pivot_longer(cols = all_of(years),
                            names_to = "year",
                            values_to = "N") |>
               mutate(N = N * selectivity_values[as.character(age)]) |>     # apply selectivity to obtain surveyed pop
               summarise(N = sum(N), .by = "year") |> # calculate the sum of N across ages
               mutate(rel_N = N/mean(N), # standardize the annual population by the average population size over the projection
                      year = as.integer(year),
                      scenario = "True")
) |>
  map_dfr(~pluck(.), .id = "pop")


trueN_stat <- trueN |>
  group_by(year) |>
  summarise(
    mean_rel_N = mean(rel_N),
    sd_rel_N = sd(rel_N),
    n_pops = n(),
    se_rel_N = sd_rel_N / sqrt(n_pops),
    cv_rel_N = sd_rel_N / mean_rel_N,
    sdlog = sqrt(log(1 + cv_rel_N^2)),
    meanlog = log(mean_rel_N) - 0.5 * sdlog^2,
    ci_lower = qlnorm(0.025, meanlog = meanlog, sdlog = sdlog),
    ci_upper = qlnorm(0.975, meanlog = meanlog, sdlog = sdlog)
  )


ggplot(trueN_stat, aes(x = year, y = mean_rel_N)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line() +
  theme_bw()





## Abundance Index ####
# calculate the abundance index and relative abundance index for each of the scenarios

# extract the strata that were used to predict spatial distributions in sdmTMB
strat <- map(dist, ~unique(.$strat))

source(here("R/stratmean_fn.R"))

stratmean_sq_all <- calc_stratmean(
  surv_list = survdat_sq,
  strata_wts = strata_wts,
  survey_area = survey_area,
  scenario_name = "Status Quo",
  value_col = "n"
)

# Group by pop (and sim if needed) to compute rel_ihat within each population
ihat_sq_all <- stratmean_sq_all |>
  group_by(pop, sim) |> #Calculations are done within each population realization (pop) and survey replicate (sim)
  mutate(
    n_years = n_distinct(year), #counts the number of unique years - needed to compute cariance of the mean across years
    mean_ihat = mean(stratmu, na.rm = TRUE), #average index across years

    # variance of the mean over years - assuming years are independent
    var_mean_ihat = sum(stratvar, na.rm = TRUE) / (n_years^2),

    # covariance between each annual stratified mean and the mean over years
    cov_stratmu_mean = stratvar / n_years,

    # relative abundance index - standardizes each year relative to the average
    rel_ihat = stratmu / mean_ihat,

    # delta-method variance for rel_ihat = stratmu / mean_ihat (Lohr, 2019 Ch9)
    rel_var =
      (stratvar / (mean_ihat^2)) +
      ((stratmu^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * stratmu * cov_stratmu_mean / (mean_ihat^3)),

    rel_se = sqrt(pmax(rel_var,0)), #standard error of estimator
    rel_cv = rel_se / rel_ihat, #coefficient of variation

    rel_log_sd = sqrt(log(1 + rel_cv^2)), #convert the CV into the lognormal SD parameter
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,

    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) %>%
  ungroup()



### Precluded Survey ####
min_tows <- 3

# 1. find strata that actually have enough tows
valid_strata <- map_dfr(survdat_precl, ~as_tibble(.x)) |>
  filter(year %in% 6:15) |>
  group_by(strat) |>
  summarise(towct = n_distinct(set), .groups = "drop") |>
  filter(towct >= min_tows) |>
  pull(strat)

# 2. filter survey data
survdat_precl_filt <- map(survdat_precl, ~
                            as_tibble(.x) |>
                            filter(strat %in% valid_strata)
)

# 3. filter weights + recompute area
strata_wts_filt <- strata_wts |>
  filter(strat %in% valid_strata)

survey_area_filt <- sum(strata_wts_filt$Area_SqNm, na.rm = TRUE)


# ---- years 1:5: same domain as status quo ----
survdat_precl_y1_5 <- map(survdat_precl, ~
                            as_tibble(.x) |>
                            filter(year %in% 1:5)
)

stratmean_precl_y1_5 <- calc_stratmean(
  surv_list     = survdat_precl_y1_5,
  strata_wts    = strata_wts,
  survey_area   = survey_area,
  scenario_name = "Preclusion",
  value_col     = "n",
  years         = 1:5
)

# ---- years 6:15: filter to strata still present under preclusion ----
survdat_precl_y6_15 <- map(survdat_precl, ~
                             as_tibble(.x) |>
                             filter(year %in% 6:15)
)

stratmean_precl_y6_15 <- calc_stratmean(
  surv_list     = survdat_precl_y6_15,
  strata_wts    = strata_wts_filt,
  survey_area   = survey_area_filt,
  scenario_name = "Preclusion",
  value_col     = "n",
  years         = 6:15
)

# ---- merge ----
stratmean_precl_all <- bind_rows(stratmean_precl_y1_5, stratmean_precl_y6_15) |>
  arrange(pop, sim, year)


# Group to compute rel_ihat
ihat_precl_all <- stratmean_precl_all |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(stratmu, na.rm = TRUE),
    var_mean_ihat = sum(stratvar, na.rm = TRUE) / (n_years^2),
    cov_stratmu_mean = stratvar / n_years,
    rel_ihat = stratmu / mean_ihat,
    rel_var =
      (stratvar / (mean_ihat^2)) +
      ((stratmu^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * stratmu * cov_stratmu_mean / (mean_ihat^3)),
    rel_se = sqrt(pmax(rel_var,0)),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) %>%
  ungroup()


### Reallocation Survey ####
min_tows <- 3

# 1. find strata that actually have enough tows
valid_strata_r <- map_dfr(survdat_reall, ~as_tibble(.x)) |>
  filter(year %in% 6:15) |>
  group_by(strat) |>
  summarise(towct = n_distinct(set), .groups = "drop") |>
  filter(towct >= min_tows) |>
  pull(strat)

# 2. filter survey data
survdat_reall_filt <- map(survdat_reall, ~
                            as_tibble(.x) |>
                            filter(strat %in% valid_strata_r)
)

# 3. filter weights + recompute area
strata_wts_filt_r <- strata_wts |>
  filter(strat %in% valid_strata_r)

survey_area_filt_r <- sum(strata_wts_filt_r$Area_SqNm, na.rm = TRUE)


# ---- years 1:5: same domain as status quo ----
survdat_reall_y1_5 <- map(survdat_reall, ~
                            as_tibble(.x) |>
                            filter(year %in% 1:5)
)

stratmean_reall_y1_5 <- calc_stratmean(
  surv_list     = survdat_reall_y1_5,
  strata_wts    = strata_wts,
  survey_area   = survey_area,
  scenario_name = "Reallocation",
  value_col     = "n",
  years         = 1:5
)

# ---- years 6:15: filter to strata still present under preclusion ----
survdat_reall_y6_15 <- map(survdat_reall, ~
                             as_tibble(.x) |>
                             filter(year %in% 6:15)
)

stratmean_reall_y6_15 <- calc_stratmean(
  surv_list     = survdat_reall_y6_15,
  strata_wts    = strata_wts_filt_r,
  survey_area   = survey_area_filt_r,
  scenario_name = "Reallocation",
  value_col     = "n",
  years         = 6:15
)

# ---- merge ----
stratmean_reall_all <- bind_rows(stratmean_reall_y1_5, stratmean_reall_y6_15) |>
  arrange(pop, sim, year)


# Group to compute rel_ihat
ihat_reall_all <- stratmean_reall_all |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(stratmu, na.rm = TRUE),
    var_mean_ihat = sum(stratvar, na.rm = TRUE) / (n_years^2),
    cov_stratmu_mean = stratvar / n_years,
    rel_ihat = stratmu / mean_ihat,
    rel_var =
      (stratvar / (mean_ihat^2)) +
      ((stratmu^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * stratmu * cov_stratmu_mean / (mean_ihat^3)),
    rel_se = sqrt(pmax(rel_var,0)),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) %>%
  ungroup()


### Bind Indices ####
# bind all three scenario dataframes for efficient plotting and statistic calculation
indices <- bind_rows(ihat_sq_all,ihat_precl_all,ihat_reall_all)


## SAVE THE DATA ####
saveRDS(trueN, here(surv.prods, str_c(species, season, "TrueN.rds", sep = "_")))
saveRDS(indices, here(surv.prods, str_c(species, season, "all-ihat-surveys.rds", sep = "_")))

saveRDS(ihat_sq_all, here(surv.prods, str_c(species, season, "100pops-25sims-sq_rel-ihat.rds", sep = "_")))
saveRDS(ihat_precl_all, here(surv.prods, str_c(species, season, "100pops-25sims-precl_rel-ihat.rds", sep = "_")))
saveRDS(ihat_reall_all, here(surv.prods, str_c(species, season, "100pops-25sims-reall_rel-ihat.rds", sep = "_")))



#Read ratio est and model based
trueN <- readRDS(here(surv.prods, str_c(species, season, "TrueN.rds", sep="_")))
indices <- readRDS(here(surv.prods, str_c(species, season, "all-ihat-surveys.rds", sep="_")))
ihat_ratioest_all <- readRDS(here(surv.prods, "ratio_est", "summerflounder", "summerflounder_fall_ratio_estimator_ihat.rds"))

ihat_model_all <- readRDS(here(surv.prods, "fit_out", "sumflounder", "summerflounder_fall_model_based_ihat.rds"))
ihat_model_wind_all <- readRDS(here(surv.prods, "fit_out", "sumflounder", "summerflounder_fall_model_based_wind_ihat.rds"))


ihat_model_all2 <- ihat_model_all |>
  rename(year = YEAR) |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(est, na.rm = TRUE),
    cv = sqrt(exp(se^2) - 1),
    var_mean_ihat = sum(est^2 * cv^2, na.rm = TRUE) / (n_years^2), # cprrected using se; sqrt(exp(se^2) - 1) = (est^2 * cv^2); this is the model-based analog of var_mean_ihat
    cov_est_mean =  (est^2 * cv^2) / n_years,  # covariance between annual estimate and the across-year mean
    rel_ihat = est / mean_ihat,
    rel_var =
      ((est^2 * cv^2) / (mean_ihat^2)) +
      ((est^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * est * cov_est_mean / (mean_ihat^3)),
    rel_var = pmax(rel_var, 0),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) |>
  ungroup()


ihat_model_final <- ihat_model_all2 |>
  mutate(stratmu = NA_real_,
         stratvar = NA_real_,
         scenario = type   # rename "type" to "scenario"
  ) |>
  rename(cov_stratmu_mean = cov_est_mean) |>
  select(pop, sim, year, stratmu, stratvar, cv, scenario, n_years, mean_ihat, var_mean_ihat, cov_stratmu_mean,
         rel_ihat, rel_var, rel_se, rel_cv, rel_log_mean, rel_log_sd, rel_ci_lower, rel_ci_upper)



ihat_model_wind_all2 <- ihat_model_wind_all |>
  rename(year = YEAR) |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(est, na.rm = TRUE),
    cv = sqrt(exp(se^2) - 1),
    var_mean_ihat = sum(est^2 * cv^2, na.rm = TRUE) / (n_years^2), # cprrected using se; sqrt(exp(se^2) - 1) = (est^2 * cv^2); this is the model-based analog of var_mean_ihat
    cov_est_mean =  (est^2 * cv^2) / n_years,  # covariance between annual estimate and the across-year mean
    rel_ihat = est / mean_ihat,
    rel_var =
      ((est^2 * cv^2) / (mean_ihat^2)) +
      ((est^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * est * cov_est_mean / (mean_ihat^3)),
    rel_var = pmax(rel_var, 0),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) |>
  ungroup()


ihat_model_wind_final <- ihat_model_wind_all2 |>
  mutate(stratmu = NA_real_,
         stratvar = NA_real_,
         scenario = type   # rename "type" to "scenario"
  ) |>
  rename(cov_stratmu_mean = cov_est_mean) |>
  select(pop, sim, year, stratmu, stratvar, cv, scenario, n_years, mean_ihat, var_mean_ihat, cov_stratmu_mean,
         rel_ihat, rel_var, rel_se, rel_cv, rel_log_mean, rel_log_sd, rel_ci_lower, rel_ci_upper)



ihat_model_final <- ihat_model_final |> mutate(scenario = "Model based")
ihat_model_wind_final <- ihat_model_wind_final |>  mutate(scenario = "Model based wind")



indices2 <- bind_rows(indices, ihat_ratioest_all)
indices_models <- bind_rows(ihat_model_final, ihat_model_wind_final)

indices_models <- indices_models |>
  select(names(indices2))


indices3 <- bind_rows(indices2, indices_models)

indices3 |> filter(sim == 1, pop ==1, year==6)

saveRDS(indices3, here(surv.prods, str_c(species, season, "indices.rds", sep = "_")))




