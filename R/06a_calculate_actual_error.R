### created: 01/18/2024
### updated:

# 06a - CALCULATE ACTUAL ERROR ####


## Objective ####
#

# Outputs:

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
# source(here("R", "sim_stratmean_fn.R"))
theme_set(theme_bw())


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


## CALCULATE RELATIVE AND ABSOLUTE ERRORS ####
errors <- indices |>
   left_join(trueN, by = c("sim", "year")) |>
   mutate(rel_err = (rel_ihat-rel_N)/rel_ihat,
          abs_rel_err = abs(rel_err)) |>
  dplyr::select(!scenario.y) |>
  rename(scenario = scenario.x)


## PLOTS ####
# relative error plot
ggplot(errors) +
  geom_boxplot(aes(x = as.factor(year), y = rel_err, color = scenario)) +
  # ylim(0, NA) +
  labs(x = "Year", y = "Relative error", title = str_c("Distribution of relative errors for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "50RelErrBoxPlot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

# absolute relative error plot
ggplot(errors) +
  geom_boxplot(aes(x = as.factor(year), y = abs_rel_err, color = scenario)) +
  ylim(0, NA) +
  labs(x = "Year", y = "Absolute relative error", title = str_c("Distribution of absolute relative errors for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "50AbsRelErrBoxPlot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)


## SAVE THE DATA ####
saveRDS(errors, here(perform.metrics, str_c(species, season, "all-50rel-error.rds", sep = "_")))

