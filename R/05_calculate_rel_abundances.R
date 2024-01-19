### created: 01/18/2024
### updated:

# 05 - CALCULATE RELATIVE ABUNDANCES ####


## Objective ####
# For a given species, distribution, and survey, calulate the relative true abundnace and the relative abundance index over time.

# Outputs:

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
source(here("R", "sim_stratmean_fn.R"))
theme_set(theme_bw())


### DATA SET UP ####
# data locations
# sseep.analysis <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis"
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

### name of species to be simulated
species <- "sumflounder"

### season to be simulated
season <- "fall"


### LOAD DATA ####
# simulated abundance and distributions created here("R", "03_append_distributions.R")
pop <- readRDS(here(dist.dat, str_c(species, "_abund-dist.rds", sep = "")))

# simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
survdat_sq <- readRDS(here(survdat, str_c(species, season, "sq-survdat.rds", sep = "_")))

# simulated precluded survey data created here("R", "04b_simulate_status-quo-survey.R")
survdat_precl <- readRDS(here(survdat, str_c(species, season, "precl-survdat.rds", sep = "_")))

# simulated reallocated survey data created here("R", "04c_simulate_status-quo-survey.R")
survdat_reall <- readRDS(here(survdat, str_c(species, season, "reall-survdat.rds", sep = "_")))

# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# find the total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))


## True Abundance ####
# calculate the relative true abundance from the simulated population and distribution
trueN <- map(pop, ~as_tibble(.$pop$N) |>
                mutate(age = ages) |>
                pivot_longer(cols = all_of(years),
                             names_to = "year",
                             values_to = "N") |>
                summarise(N = sum(N), .by = "year") |> # calculate the sum of N across ages
                mutate(rel_N = N/mean(N), # standardize the annual population by the average population size over the projection
                       year = as.integer(year),
                       scenario = "True")
              ) |>
  map_dfr(~pluck(.), .id = "sim")

# calculate the median relative abundance value and the upper and lower confidence intervals across simulations
true_med <- trueN |>
  group_by(year, scenario) |>
  summarise(med = median(rel_N),
            lower = quantile(rel_N, 0.025),
            upper = quantile(rel_N, 0.975)) |>
  mutate(type = "Relative True Abundance")



## Abundance Index ####
# calculate the abundance index and relative abundance index for each of the scenarios

### Status Quo ####
ihat_sq <- map(survdat_sq, ~as_tibble(.$setdet) |> sim_stratmean( strata_wts = strata_wts, survey_area = survey_area) |>
              mutate(rel_ihat = stratmu/mean(stratmu),
                     scenario = "Status Quo")) |>
  map_dfr(~pluck(.), .id = "sim")

### Precluded Survey ####
ihat_precl <- map(survdat_precl, ~as_tibble(.) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
                 mutate(rel_ihat = stratmu/mean(stratmu),
                        scenario = "Preclusion")) |>
  map_dfr(~pluck(.), .id = "sim")

### Reallocated Survey ####
ihat_reall <- map(survdat_reall, ~as_tibble(.x) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
                 mutate(rel_ihat = stratmu/mean(stratmu),
                        scenario = "Reallocation")) |>
  map_dfr(~pluck(.), .id = "sim")

### Bind Indices ####
# bind all three scenario dataframes for efficient plotting and statistic calculation
indices <- bind_rows(ihat_sq, ihat_precl, ihat_reall)

# calculate the median relative abundance indices and the upper and lower confidence intervals across simulations and scenarios
indices_med <- indices |>
  group_by(year, scenario) |>
  summarise(med = median(rel_ihat),
            lower = quantile(rel_ihat, 0.025),
            upper = quantile(rel_ihat, 0.975)) |>
  mutate(type = "Estimated Relative Abundance Index")


## Plots ####
# bind the relative data frames and plot to compare scenarios to relative true abudnance
bind_rows(true_med, indices_med) |>
ggplot() +
  aes(x = year) +
  geom_line(aes(y = med, color = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = scenario), alpha = 0.25) +
  #geom_line(aes(y = rel_N), color = "red") +
  labs(x = "Year", y = "Relative Abundance ", color = "Scenario", fill = "Scenario", title = str_c("Comparison of median relative abundance across scenarios for", season, species, "populations", sep = " ")) +
  # scale_color_manual(name = "Type", values = c("darkblue", "red")) +
  # scale_fill_manual(name = "Type", values = c("darkblue", "red")) +
  theme(legend.position = "bottom")

ggsave(str_c(species, season, "MedRelAbundPlot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 10, height = 6)


## SAVE THE DATA ####
saveRDS(trueN, here(surv.prods, str_c(species, season, "rel-TrueN.rds", sep = "_")))
saveRDS(indices, here(surv.prods, str_c(species, season, "all-ihat.rds", sep = "_")))
saveRDS(ihat_sq, here(surv.prods, str_c(species, season, "sq_rel-ihat.rds", sep = "_")))
saveRDS(ihat_precl, here(surv.prods, str_c(species, season, "precl_rel-ihat.rds", sep = "_")))
saveRDS(ihat_reall, here(surv.prods, str_c(species, season, "reall_rel-ihat.rds", sep = "_")))

