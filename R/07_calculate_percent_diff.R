### created: 02/05/2024
### updated: 02/06/2024

# 07 - CALCULATE PERCENT DIFFERENCE FOR SCENARIOS ####


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

### ages simulated
ages <- 0:7

### years projected
years <- 1:5

### number of simulations
nsims <- 1:2

### number of simulations
nsurveys <- 20

### LOAD DATA ####
# relative abunance indices for each survey effort scenario created here("R", "05_calculate_rel_abundance.R")
### STATUS QUO
# ihat_sq <- readRDS(here(surv.prods, str_c(species, season, "sq_rel-ihat.rds", sep = "_")))
ihat_sq1 <- readRDS(here(surv.prods, str_c(species, season, "1pop-20sq_rel-ihat.rds", sep = "_")))


### PRECLUSION
# ihat_precl <- readRDS(here(surv.prods, str_c(species, season, "precl_rel-ihat.rds", sep = "_")))
ihat_precl1 <- readRDS(here(surv.prods, str_c(species, season, "1pop-20precl_rel-ihat.rds", sep = "_")))

### REALLOCATION
# ihat_reall <- readRDS(here(surv.prods, str_c(species, season, "reall_rel-ihat.rds", sep = "_")))
ihat_reall1 <- readRDS(here(surv.prods, str_c(species, season, "1pop-20reall_rel-ihat.rds", sep = "_")))

## CALCULATE ABSOLUTE PERCENT DIFFERENCES ####
### Status Quo vs Preclusion ####
sq_precl_diff <- left_join(ihat_sq1, ihat_precl1, by = c("sim", "year")) |>
  janitor::clean_names() |>
  group_by(sim, year) |>
  mutate(perc_diff_mu =( abs(rel_ihat_x - rel_ihat_y) / rel_ihat_x)*100,
         perc_diff_cv =( abs(cv_x - cv_y) / cv_x)*100) |>
  select(sim, year, perc_diff_mu, perc_diff_cv) |>
  mutate(type = "SQ-Preclusion")



### Status Quo vs Reallocation ####
sq_reall_diff <- left_join(ihat_sq1, ihat_reall1, by = c("sim", "year")) |>
  janitor::clean_names() |>
  group_by(sim, year) |>
  mutate(perc_diff_mu =( abs(rel_ihat_x - rel_ihat_y) / rel_ihat_x)*100,
         perc_diff_cv =( abs(cv_x - cv_y) / cv_x)*100) |>
  select(sim, year, perc_diff_mu, perc_diff_cv) |>
  mutate(type = "SQ-Reallocation")


## BIND ####
diffs <- bind_rows(sq_precl_diff, sq_reall_diff)

## PLOT DISTRIBUTION ####
ggplot(diffs) +
  geom_boxplot(aes(x = as.factor(year), y = perc_diff_mu, color = type)) +
  labs(x = "Year", y = "Absolute percent difference", subtitle = "Distribution of the absolute percent difference between\nrelative abundance indices when compared to the status quo index")

ggsave(str_c(species, season, "perc-diff-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)
