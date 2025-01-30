### created: 01/22/2024
### updated: 02/05/2024

# 06c - CALCULATE ERROR IN TREND OVER TIME ####


## Objective ####
#

# Outputs:

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
library(broom)
# source(here("R", "sim_stratmean_fn.R"))


### DATA SET UP ####
# data locations
# sseep.analysis <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis"
# survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")
plots <- here("outputs", "plots")

### name of species to be simulated
species <- "scup"

### season to be simulated
season <- "fall"

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### number of simulations
nsims <- 1:2

### number of simulations
nsurveys <- 25

### LOAD DATA ####
# relative true abundance created here("R", "05_calculate_rel_abundance.R")
trueN <- readRDS(here(surv.prods, str_c(species, season, "rel-TrueN.rds", sep = "_")))

# relative abunance indices across scenarios created here("R", "05_calculate_rel_abundance.R")
# indices <- readRDS(here(surv.prods, str_c(species, season, "all-ihat.rds", sep = "_")))
indices25 <- readRDS(here(surv.prods, str_c(species, season, "all-ihat_1pop-25survs.rds", sep = "_")))


## ERROR IN TREND OVER TIME ####
# relative true abundance trends
trueN_lm <- trueN |>
  group_by(sim, scenario) |>
  nest() |>
  mutate(mods = map(data, ~lm(rel_N ~ year, data = .x)),
         coefs = map(mods, broom::tidy, conf.int = TRUE),
         slope = map(coefs, ~filter(.x, term == "year"))) |>
  select(sim, slope, scenario) |>
  unnest(cols = slope)

saveRDS(trueN_lm, here(surv.prods, str_c(species, season, "TrueRelAbundTrends.rds", sep = "_")))

indices25_lm <- indices25 |>
  group_by(sim, scenario) |>
  nest() |>
  mutate(mods = map(data, ~lm(rel_ihat ~ year, data = .x)),
         coefs = map(mods, broom::tidy, conf.int = TRUE),
         slope = map(coefs, ~filter(.x, term == "year"))) |>
  select(sim, slope, scenario) |>
  unnest(cols = slope)

saveRDS(indices25_lm, here(surv.prods, str_c(species, season, "RelIndexTrends.rds", sep = "_")))


slope_errors <- indices25_lm |>
  select(sim, term, estimate) |>
  # left_join(trueN_lm[1:3], by = "sim") |>
  bind_cols(trueN_lm[1,1:3]) |>
  janitor::clean_names() |>
  mutate(err = estimate_4-estimate_7,
         rel_err = (err)/estimate_7,
         abs_rel_err = abs(rel_err))

saveRDS(slope_errors, here(surv.prods, str_c(species, season, "SlopeErrors.rds", sep = "_")))
#|>
  # group_by(scenario) |>
  # summarise(med_rel_err = median(rel_err),
  #           med_abs_err = median(abs_rel_err),
  #           lower_rel = quantile(rel_err, 0.025),
  #           upper_rel = quantile(rel_err, 0.975),
  #           lower_abs = quantile(abs_rel_err, 0.025),
  #           upper_abs = quantile(abs_rel_err, 0.975)) #|>
  #   ggplot() +
  # aes(x = as.factor(scenario), y = med_rel_err) +
  # geom_point(position = position_dodge2(width = 0.2)) +
  # geom_errorbar(aes(ymin=lower_rel, ymax=upper_rel), width= 0.2, position = position_dodge2(width = 0.2))


## PLOTS ####

# ggplot(indices_lm) +
#   geom_boxplot(aes(x = as.factor(scenario), y = estimate, color = scenario)) +
#   # ylim(0, NA) +
#   labs(x = "Year", y = "Linear regression slope estimates", title = str_c("Distribution of linear regression slope estimates for", season, species, "survey", sep = " ")) +
#   theme(legend.position = "bottom")

ggplot() +
  geom_boxplot(data = indices25_lm, aes(x = fct_inorder(scenario), y = estimate)) +
  geom_hline(yintercept = trueN_lm$estimate[1], lty = 5, color = "red", size = 0.75, show.legend = TRUE) +
  annotate("text", y = 0.25, x = 2, label = str_c("Relative true abunandance slope:", round(trueN_lm$estimate[1], 2), sep = " "), color = "red") +
  labs(x = "Scenario", y = "Linear regression slope estimate", subtitle = str_c("Distribution of linear regression slope estimates for", season, species, "survey indices over time", sep = " "))

ggsave(str_c(species, season, "lm-ests-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)


ggplot(slope_errors) +
  geom_boxplot(aes(x = as.factor(scenario), y = err, color = fct_inorder(scenario))) +
  # ylim(0, NA) +
  labs(x = "Year", y = "Linear regression slope error", title = str_c("Distribution of errors of linear regression slope estimates for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "lm-errors-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

ggplot(slope_errors) +
  geom_boxplot(aes(x = as.factor(scenario), y = rel_err, color = fct_inorder(scenario))) +
  # ylim(0, NA) +
  labs(x = "Year", y = "Linear regression slope relative errors", title = str_c("Distribution of relative errors of linear regression slope estimates for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "lm-rel-errors-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

ggplot(slope_errors) +
  geom_boxplot(aes(x = as.factor(scenario), y = abs_rel_err, color = fct_inorder(scenario))) +
  # ylim(0, NA) +
  labs(x = "Year", y = "Linear regression slope absolute relative errors", title = str_c("Distribution of absolute relative errors of linear regression slope estimates for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "lm-abs-rel-errors-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 9, height = 6)

