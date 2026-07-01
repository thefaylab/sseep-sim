# ============================================================
# Supplemental correction factor estimator
# Random stratified BTS mean + supplemental sampling correction
# ============================================================

library(tidyverse)
library(here)
library(sdmTMB)
library(kableExtra)

source(here("R", "stratmean_fn.R"))

# ----------------------------
# DATA SET UP
# ----------------------------

sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

species <- "scup"
season <- "fall"
ages <- 0:7
years <- 1:15
nsims <- 1:100
ids <- sprintf("%03d", nsims)

min_tows <- 3

# ----------------------------
# LOAD DATA
# Loads three survey scenarios
# ----------------------------

survdat_sq <- map(ids, function(id) {
  x <- readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, id)))
  out <- x$setdet
  rm(x); gc()
  out
})

survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))

survdat_supp <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_supp+precl_survey.rds", species, season, .x))))


dist <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))

strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

survey_area <- sum(strata_wts$Area_SqNm, na.rm = TRUE)

strat <- map(dist, ~unique(.x$strat))

# ----------------------------
# STATUS QUO INDEX
# ----------------------------

stratmean_sq_all <- calc_stratmean(
  surv_list = survdat_sq,
  strata_wts = strata_wts,
  survey_area = survey_area,
  scenario_name = "Status Quo",
  value_col = "n"
)

# ----------------------------
# BTS OUTSIDE-WIND MEANS
# years 6-15, random stratified BTS when preclusion starts
# ----------------------------

bts_out_mu <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  as_tibble(surv) |>
    filter(strat %in% unlist(strat), year %in% 6:15, AREA_CODE != 1) |>
    group_by(pop = pop_num, sim, year, strat) |>
    summarise(
      towct_out = n_distinct(set),
      mu_out = sum(n, na.rm = TRUE) / towct_out,
      var_out = ifelse(towct_out < 2, 0, var(n, na.rm = TRUE)),
      .groups = "drop"
    )
})

# ----------------------------
# SUPPLEMENTAL INSIDE-WIND MEANS
# years 6-15, supplemental survey inside wind areas
# ----------------------------

supp_in_mu <- map2_dfr(survdat_supp, seq_along(survdat_supp), function(surv, pop_num) {
  as_tibble(surv) |>
    filter(strat %in% unlist(strat), year %in% 6:15, AREA_CODE == 1) |>
    group_by(pop = pop_num, sim, year, strat) |>
    summarise(
      towct_in = n_distinct(set),
      mu_in = sum(n, na.rm = TRUE) / towct_in,
      var_in = ifelse(towct_in < 2, 0, var(n, na.rm = TRUE)),
      .groups = "drop"
    )
})

# ----------------------------
# SUPPLEMENTAL CORRECTION FACTOR
# ratio_supp = inside-wind mean / outside-wind mean
# ----------------------------

supp_ratio_tbl <- bts_out_mu |>
  left_join(supp_in_mu, by = c("pop", "sim", "year", "strat")) |>
  mutate(
    ratio_supp = if_else(
      towct_out >= min_tows & towct_in >= min_tows &
        !is.na(mu_out) & !is.na(mu_in) & mu_out > 0,
      mu_in / mu_out,
      NA_real_
    )
  )

# ----------------------------
# APPLY CORRECTION FACTOR
# corrected mean = BTS outside-wind mean * supplemental correction factor
# ----------------------------

corrected_supp_post <- supp_ratio_tbl |>
  filter(!is.na(ratio_supp)) |>
  left_join(strata_wts, by = "strat") |>
  mutate(
    W = Area_SqNm / survey_area,
    mu_corrected = mu_out * ratio_supp,
    var_corrected = ratio_supp^2 * var_out,
    wt_mu = W * mu_corrected,
    wt_var = ifelse(towct_out < 2, 0, (W^2) * var_corrected / towct_out)
  ) |>
  group_by(pop, sim, year) |>
  summarise(
    stratmu = sum(wt_mu, na.rm = TRUE),
    stratvar = sum(wt_var, na.rm = TRUE),
    cv = sqrt(stratvar) / stratmu,
    .groups = "drop"
  ) |>
  mutate(scenario = "Supplemental correction factor")

# ----------------------------
# YEARS 1-5 COME FROM STATUS QUO
# no wind preclusion yet
# ----------------------------

corrected_supp_pre <- stratmean_sq_all |>
  filter(year %in% 1:5) |>
  select(pop, sim, year, stratmu, stratvar, cv, scenario) |>
  mutate(scenario = "Supplemental correction factor")

corrected_supp_all <- bind_rows(corrected_supp_pre, corrected_supp_post) |>
  arrange(pop, sim, year)

# ----------------------------
# RELATIVE INDEX + UNCERTAINTY
# ----------------------------

supp_corr_ihat <- corrected_supp_all |>
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
    rel_se = sqrt(pmax(rel_var, 0)),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) |>
  ungroup()

# ----------------------------
# SAVE
# ----------------------------

supp_corr.dir <- here("data", "rds", "surv-prods", "supp_corr", species)
dir.create(supp_corr.dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  supp_corr_ihat,
  here(supp_corr.dir, str_c(species, season, "supplemental_correction_ihat.rds", sep = "_"))
)
