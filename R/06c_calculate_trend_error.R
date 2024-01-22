### created: 01/22/2024
### updated:

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
surv.prod <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")
plots <- here("outputs", "plots")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"

### LOAD DATA ####
# relative true abundance created here("R", "05_calculate_rel_abundance.R")
trueN <- readRDS(here(surv.prods, str_c(species, season, "rel-TrueN.rds", sep = "_")))

# relative abunance indices across scenarios created here("R", "05_calculate_rel_abundance.R")
indices <- readRDS(here(surv.prods, str_c(species, season, "all-ihat.rds", sep = "_")))


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


indices_lm <- indices |>
  group_by(sim, scenario) |>
  nest() |>
  mutate(mods = map(data, ~lm(rel_ihat ~ year, data = .x)),
         coefs = map(mods, broom::tidy, conf.int = TRUE),
         slope = map(coefs, ~filter(.x, term == "year"))) |>
  select(sim, slope, scenario) |>
  unnest(cols = slope)


indices_lm |>
  select(sim, term, estimate) |>
  left_join(trueN_lm[1:3], by = "sim") |>
  janitor::clean_names() |>
  mutate(rel_err = (estimate_x-estimate_y)/estimate_x,
         abs_rel_err = abs(rel_err)) |>
  group_by(scenario) |>
  summarise(med_rel_err = median(rel_err),
            med_abs_err = median(abs_rel_err),
            lower_rel = quantile(rel_err, 0.025),
            upper_rel = quantile(rel_err, 0.975),
            lower_abs = quantile(abs_rel_err, 0.025),
            upper_abs = quantile(abs_rel_err, 0.975)) #|>
  #   ggplot() +
  # aes(x = as.factor(scenario), y = med_rel_err) +
  # geom_point(position = position_dodge2(width = 0.2)) +
  # geom_errorbar(aes(ymin=lower_rel, ymax=upper_rel), width= 0.2, position = position_dodge2(width = 0.2))




ggplot(indices_lm) +
  geom_boxplot(aes(x = as.factor(scenario), y = estimate, color = scenario)) +
  # ylim(0, NA) +
  labs(x = "Year", y = "CV", title = str_c("Distribution of CVs for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "50CV-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

