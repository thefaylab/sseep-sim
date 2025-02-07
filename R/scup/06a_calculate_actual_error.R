### created: 01/18/2024
### updated: 02/07/2025

# 06a - CALCULATE ACTUAL ERROR ####


## Objective ####
# For a given species and iteration, calculate the relative and absolute error to compare an abundance index to the true abundance

# Outputs:
#  a dataframe with values for each simulation and year pertaining to the relative error and absolute relative error of a given survey and its abundance index compared to the relative true abundance
# the relative errors and absolute errors for a given survey scenario plotted across time

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
years <- 1:15

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

## CALCULATE RELATIVE AND ABSOLUTE ERRORS ####
errors25 <- indices25 |> mutate(sim = as.character(sim)) |>
   left_join(trueN, by = "year") |>
   mutate(rel_err = (rel_ihat-rel_N)/rel_N,
          abs_rel_err = abs(rel_err)) |>
  dplyr::select(!scenario.y) |> #delete TRUE name column only.
  rename(scenario = scenario.x) #scenarios are relative to TRUE rel N


## PLOTS ####
# relative error plot
ggplot(errors25) +
  geom_boxplot(aes(x = as.factor(year), y = rel_err, color = fct_inorder(scenario))) +
  # ylim(0, NA) +
  labs(x = "Year", y = "Relative error", title = str_c("Distribution of relative errors for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "RelErrBoxPlot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

# absolute relative error plot
ggplot(errors25) +
  geom_boxplot(aes(x = as.factor(year), y = abs_rel_err, color = fct_inorder(scenario))) +
  ylim(0, NA) +
  labs(x = "Year", y = "Absolute relative error", title = str_c("Distribution of absolute relative errors for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "AbsRelErrBoxPlot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)


## SAVE THE DATA ####
saveRDS(errors25, here(perform.metrics, str_c(species, season, "1pop_25all-rel-error.rds", sep = "_")))

