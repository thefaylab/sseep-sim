### created: 01/18/2024
### updated:

# 06b - CALCULATE ESTIMATION ERROR ####


## Objective ####
# For a given species, distribution, and survey, calulate

# Outputs:

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
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

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### number of simulations
nsims <- 1:2

### number of simulations
nsurveys <- 20

### LOAD DATA ####
# relative and absolute errors of simulated abundance created here("R", "06a_calculate_actual_abundance.R")
errors <- readRDS(here(perform.metrics, str_c(species, season, "all-rel-error.rds", sep = "_")))


## ESTIMATION ERROR ####
ggplot(errors) +
  geom_boxplot(aes(x = as.factor(year), y = cv, color = fct_inorder(scenario))) +
  ylim(0, NA) +
  labs(x = "Year", y = "CV", title = str_c("Distribution of CVs for", season, species, "survey", sep = " ")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "CV-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

