# OBJECTIVE

#For a given species, calculate a ratio estimator as a potential mitigation approach to adjust for the loss of sampling coverage caused by offshore wind energy development.

#Outputs: abundance indices adjusted with the ratio for each population and simulation of the projection.


# PACKAGES

library(tidyverse)
library(here)
library(sdmTMB)
source(here("R", "sim_stratmean_fn.R"))
library(kableExtra)

# DATA SET UP
## Directories
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

## Parameters
species <- "scup"
season  <- "fall"
ages      <- 0:7
years     <- 1:15
nsims   <- 1:100
ids     <- sprintf("%03d", nsims)

## Data
#survdat_sq <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, .x))))
survdat_sq <- map(ids, function(id) {
  x <- readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds",
                                     species, season, id)))
  out <- x$setdet #loas only setdet data
  rm(x); gc()
  out
})

survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))
dist <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))


# Area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# Total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))

# ABUDNANCE INDEX
#calculate the abundance index and relative abundance index for each of the scenarios

## Status Quo Survey

#extract the strata that were used to predict spatial distributions in sdmTMB
strat <- map(dist, ~unique(.$strat))

source(here("R/stratmu_fn.R"))

ihat_sq_all <- calc_ihat(
  surv_list = survdat_sq,
  strata_wts = strata_wts,
  survey_area = survey_area,
  scenario_name = "Status Quo",
  value_col = "n"
)

# Group by pop (and sim if needed) to compute rel_ihat within each population
ihat_sq_all <- ihat_sq_all %>%
  group_by(pop, sim) %>%
  mutate(
    mean_ihat = mean(stratmu),
    rel_ihat = stratmu / mean_ihat,
    rel_var = stratvar / (mean_ihat^2),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_se = sqrt(log(1 + rel_cv^2)),
    rel_ci_lower = exp(log(rel_ihat) - 1.96 * rel_log_se),
    rel_ci_upper = exp(log(rel_ihat) + 1.96 * rel_log_se)
  ) %>%
  ungroup()



## Status quo stratum means (baseline)
sq_mu_base <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    filter(strat %in% unlist(strat), year %in% 1:5) |>
    group_by(pop = pop_num, sim, strat) |>
    summarise(
      towct = n_distinct(set),
      mu_sq = sum(n) / towct,
      .groups = "drop"
    )
})


## Preclusion Survey
ihat_precl_all <- calc_ihat(
  surv_list = survdat_precl,
  strata_wts = strata_wts,
  survey_area = survey_area,
  scenario_name = "Preclusion",
  value_col = "n"
)

# Group by pop (and sim if needed) to compute rel_ihat within each population
ihat_precl_all <- ihat_precl_all %>%
  group_by(pop, sim) %>%
  mutate(
    mean_ihat = mean(stratmu),
    rel_ihat = stratmu / mean_ihat,
    rel_var = stratvar / (mean_ihat^2),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_se = sqrt(log(1 + rel_cv^2)),
    rel_ci_lower = exp(log(rel_ihat) - 1.96 * rel_log_se),
    rel_ci_upper = exp(log(rel_ihat) + 1.96 * rel_log_se)
  ) %>%
  ungroup()



## Preclusion stratum means for baseline
precl_mu_y <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    filter(strat %in% unlist(strat), year %in% 1:5, AREA_CODE != 1) |>
    group_by(pop = pop_num, sim, year, strat) |>
    summarise(
      towct = n_distinct(set),
      mu_precl = sum(n) / towct,
      .groups = "drop"
    )
})


# Ratio calculation (y_tot/y_preclusion)
# ratio_tbl <- sq_mu |>
#   left_join(precl_mu, by = c("pop", "sim", "strat"))
#
# ratio_tbl <- ratio_tbl |>
#   mutate(ratio = ifelse(is.na(mu_precl) | mu_precl == 0, 1, mu_tot / mu_precl))
#



ratio_tbl <- sq_mu_y |>
  group_by(pop, sim, year, strat) |>
  summarise(mean_sq = mean(mu_sq, na.rm = TRUE), .groups = "drop") |>
  left_join(
    precl_mu_y |>
      group_by(pop, sim, year, strat) |>
      summarise(mean_precl = mean(mu_precl, na.rm = TRUE), .groups = "drop"),
    by = c("pop","sim","year","strat")
  ) |>
  mutate(ratio = ifelse(is.na(mean_precl) | mean_precl == 0, 1, mean_sq / mean_precl))



# Calculate expanded index

survdat_precl_ratio <- survdat_precl

for (i in seq_along(survdat_precl_ratio)) {
  survdat_precl_ratio[[i]] <- survdat_precl_ratio[[i]] |>
    as_tibble() |>
    mutate(pop = i) |>
    left_join(
      ratio_tbl |> select(pop, sim, year, strat,ratio),
      by = c("pop", "sim", "year","strat")
    ) |>
    mutate(
      ratio = ifelse(is.na(ratio), 1, ratio),
      n_ratio = n * ratio
    )
}


ihat_precl_exp_all <- calc_ihat(
  surv_list = survdat_precl_ratio,
  strata_wts = strata_wts,
  survey_area = survey_area,
  scenario_name = "Preclusion_expanded",
  value_col = "n_ratio",
  years = 6:15
)

ihat_precl_exp_all <- ihat_precl_exp_all %>%
  group_by(pop, sim) %>%
  mutate(
    mean_ihat = mean(stratmu),
    rel_ihat = stratmu / mean_ihat,
    rel_var = stratvar / (mean_ihat^2),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_se = sqrt(log(1 + rel_cv^2)),
    rel_ci_lower = exp(log(rel_ihat) - 1.96 * rel_log_se),
    rel_ci_upper = exp(log(rel_ihat) + 1.96 * rel_log_se)
  ) %>%
  ungroup()







#Stratum-level ratio table
sq_strat <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    mutate(pop = pop_num) |>
    group_by(pop, sim, year, strat) |>
    summarise(
      towct_sq = length(unique(set)),
      mu_sq = sum(n) / towct_sq,
      .groups = "drop"
    )
})

precl_strat <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    mutate(pop = pop_num) |>
    group_by(pop, sim, year, strat) |>
    summarise(
      towct_precl = length(unique(set)),
      mu_precl = sum(n) / towct_precl,
      .groups = "drop"
    )
})

ratio_tbl_strat <- sq_strat |>
  left_join(precl_strat, by = c("pop", "sim", "year", "strat")) |>
  mutate(
    ratio = case_when(
      is.na(mu_precl) ~ 1,
      mu_precl == 0 ~ 1,
      TRUE ~ mu_sq / mu_precl
    )
  )

survdat_precl_ratio <- survdat_precl

for (i in seq_along(survdat_precl_ratio)) {
  survdat_precl_ratio[[i]] <- survdat_precl_ratio[[i]] |>
    as_tibble() |>
    mutate(pop = i) |>
    left_join(
      ratio_tbl_strat |>
        select(pop, sim, year, strat, ratio),
      by = c("pop", "sim", "year", "strat")
    ) |>
    mutate(
      ratio = ifelse(is.na(ratio), 1, ratio),
      n_ratio = n * ratio
    )
}

ihat_precl_exp_all2 <- calc_ihat(
  surv_list = survdat_precl_ratio,
  strata_wts = strata_wts,
  survey_area = survey_area,
  scenario_name = "Preclusion_expanded",
  value_col = "n_ratio",
  years = 6:15
)


ihat_precl_exp_all2 <- ihat_precl_exp_all2 %>%
  group_by(pop, sim) %>%
  mutate(
    mean_ihat = mean(stratmu),
    rel_ihat = stratmu / mean_ihat,
    rel_var = stratvar / (mean_ihat^2),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_se = sqrt(log(1 + rel_cv^2)),
    rel_ci_lower = exp(log(rel_ihat) - 1.96 * rel_log_se),
    rel_ci_upper = exp(log(rel_ihat) + 1.96 * rel_log_se)
  ) %>%
  ungroup()
