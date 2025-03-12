### created: 01/18/2024
### updated: 02/07/2024

# 05 - CALCULATE RELATIVE ABUNDANCES ####


## Objective ####
# For a given species, distribution, and survey, calculate the relative true abundance and the relative abundance index.

# Outputs: Relative true abundance and abundance indices for each year and simulation of the projection standardized to the average abundance or abundance index over time

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
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

### name of species to be simulated
species <- "scup"

### season to be simulated
season <- "fall"

### number of simulations
nsims <- 1:2

### number of simulations
nsurveys <- 25

### ages simulated
ages <- 0:7

### years projected
years <- 1:15


### LOAD DATA ####
# simulated abundance and distributions created here("R", "03_append_distributions.R")
pop <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "abund-dist.rds", sep = "_")))
# pop_nw, here(dist.dat, str_c(species, season, length(nsims), "nowind-abund-dist.rds", sep = "_")))

# filtered sdmTMB distribution predictions created here("R", "03_append_distributions.R"); to be used to filter sampling spatial frame
dist <- readRDS(here(dist.dat, str_c(species, season, length(nsims), "sdm-dist-only.rds", sep = "_")))
# dist_nowind <- here(dist.dat, str_c(species, season, length(nsims), "sdm-dist_nowind-only.rds", sep = "_")))

# simulated status quo survey data created here("R", "04a_simulate_status-quo-survey.R")
survdat_sq <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat.rds", sep = "_")))
# survdat_sq_nw, here(survdat, str_c(species, season, length(nsims), "sims", "sq-surv-dat_nw.rds", sep = "_")))

# simulated precluded survey data created here("R", "04b_simulate_status-quo-survey.R")
survdat_precl <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims",nsurveys,  "precl-surv-dat.rds", sep = "_")))
# survdat_precl_nw, here(survdat, str_c(species, season, length(nsims), nsurveys, "precl-surv-dat_nw.rds", sep = "_")))

# simulated reallocated survey data created here("R", "04c_simulate_status-quo-survey.R")
survdat_reall <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", "reall", nsurveys, "survdat.rds", sep = "_")))
# survdat_reall_nw, here(survdat, str_c(species, season, length(nsims), "sims", "reall", nsurveys, "survdat_nw.rds", sep = "_")))




# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# find the total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))


#function to set selectivity
source(here("R/selectivity_fns.R"))
q = force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1)

ages <- 0:7 #evaluate the function over a range of ages
(selectivity_values <- q(ages))


## True Abundance ####
# calculate the relative true abundance from the simulated population and distribution
trueN <- map(pop, ~as_tibble(.$pop$N) |>
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

# calculate the median relative abundance value and the upper and lower confidence intervals across simulations
# true_med <- trueN |>
#   group_by(year, scenario) |>
#   summarise(med = median(rel_N),
#             lower = quantile(rel_N, 0.025),
#             upper = quantile(rel_N, 0.975)) |>
#   mutate(type = "Relative True Abundance")



## Abundance Index ####
# calculate the abundance index and relative abundance index for each of the scenarios

# extract the strata that were used to predict spatial distributions in sdmTMB
strat <- map(dist, ~unique(.$strat))

### Status Quo ####
# ihat_sq <- map(survdat_sq, ~as_tibble(.$setdet) |> filter(strat %in% strat) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
#               mutate(rel_ihat = stratmu/mean(stratmu),
#                      scenario = "Status Quo")) |>
#   map_dfr(~pluck(.), .id = "sim")

#FOR ONE REALIZATION OF THE POPULATION
# ihat_sq1 <- survdat_sq[[1]]$setdet |> as_tibble() |>
#   filter(strat %in% strat) |>
#   group_by(sim, year, strat) |>
#   summarise(towct = length(unique(set)),
#             mu = sum(n)/towct,
#             var = ifelse(towct == 1, 0,
#                          sum((n - mu)^2)/(towct - 1))) |>
#   left_join(strata_wts, by = "strat") |>
#   mutate(wt_mu = Area_SqNm * mu,
#          wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
#   ungroup() |>
#   group_by(sim, year) |>
#   summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
#             stratvar = sum(wt_var),
#             cv = sqrt(stratvar)/stratmu) |>
#   mutate(rel_ihat = stratmu/mean(stratmu),
#          scenario = "Status Quo")

#FOR ALL THE REALIZATIONS OF THE POPULATION
ihat_sq_all <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
  surv$setdet |>
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(towct = length(unique(set)),
              mu = sum(n)/towct,
              var = ifelse(towct == 1, 0,
                           sum((n - mu)^2)/(towct - 1)),
              .groups = "drop") |>
    left_join(strata_wts, by = "strat") |>
    mutate(wt_mu = Area_SqNm * mu,
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
    group_by(sim, year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area,
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu,
              .groups = "drop") |>
    mutate(scenario = "Status Quo",
           pop = pop_num)
})

# group by pop and sims to compute rel_ihat within each population. If it get computed in the previous step, all populations get combined so if I compute rel_ihat as stratmu / mean(stratmu) across all rows, it end up dividing by a different (larger) mean value.
ihat_sq_all <- ihat_sq_all %>%
  group_by(pop, sim) %>%
  mutate(rel_ihat = stratmu / mean(stratmu)) %>%
  ungroup()


### Precluded Survey ####
### #FOR ONE REALIZATION OF THE POPULATION
# ihat_precl <- map(survdat_precl, ~as_tibble(.) |> filter(strat %in% strat) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
#                  mutate(rel_ihat = stratmu/mean(stratmu),
#                         scenario = "Preclusion")) |>
#   map_dfr(~pluck(.), .id = "sim")

# ihat_precl1 <- survdat_precl[[1]] |> as_tibble() |>
#   filter(strat %in% strat) |>
#   group_by(sim, year, strat) |>
#   summarise(towct = length(unique(set)),
#             mu = sum(n)/towct,
#             var = ifelse(towct == 1, 0,
#                          sum((n - mu)^2)/(towct - 1))) |>
#   left_join(strata_wts, by = "strat") |>
#   mutate(wt_mu = Area_SqNm * mu,
#          wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
#   ungroup() |>
#   group_by(sim, year) |>
#   summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
#             stratvar = sum(wt_var),
#             cv = sqrt(stratvar)/stratmu) |>
#   mutate(rel_ihat = stratmu/mean(stratmu),
#          scenario = "Preclusion")

#FOR ALL THE REALIZATIONS OF THE POPULATION
ihat_precl_all <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  surv |> #it is already a data.table
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(towct = length(unique(set)),
              mu = sum(n)/towct,
              var = ifelse(towct == 1, 0,
                           sum((n - mu)^2)/(towct - 1)),
              .groups = "drop") |>
    left_join(strata_wts, by = "strat") |>
    mutate(wt_mu = Area_SqNm * mu,
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
    group_by(sim, year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area,
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu,
              .groups = "drop") |>
    mutate(scenario = "Preclusion",
           pop = pop_num)
})

# Now group by pop (and sim if needed) to compute rel_ihat within each population
ihat_precl_all <- ihat_precl_all %>%
  group_by(pop, sim) %>%
  mutate(rel_ihat = stratmu / mean(stratmu)) %>%
  ungroup()




### Reallocated Survey ####
# ihat_reall <- map(survdat_reall, ~as_tibble(.x) |> filter(strat %in% strat) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
#                  mutate(rel_ihat = stratmu/mean(stratmu),
#                         scenario = "Reallocation")) |>
#   map_dfr(~pluck(.), .id = "sim")
#

#FOR ONE REALIZATION OF THE POPULATION
# ihat_reall1 <- survdat_reall[[1]] |> as_tibble() |>
#   filter(strat %in% strat) |>
#   group_by(sim, year, strat) |>
#   summarise(towct = length(unique(set)),
#             mu = sum(n)/towct,
#             var = ifelse(towct == 1, 0,
#                          sum((n - mu)^2)/(towct - 1))) |>
#   left_join(strata_wts, by = "strat") |>
#   mutate(wt_mu = Area_SqNm * mu,
#          wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
#   ungroup() |>
#   group_by(sim, year) |>
#   summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
#             stratvar = sum(wt_var),
#             cv = sqrt(stratvar)/stratmu) |>
#   mutate(rel_ihat = stratmu/mean(stratmu),
#          scenario = "Reallocation")

#FOR ALL THE REALIZATIONS OF THE POPULATION
ihat_reall_all <- map2_dfr(survdat_reall, seq_along(survdat_reall), function(surv, pop_num) {
  surv |> #it is already a data.table
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(towct = length(unique(set)),
              mu = sum(n)/towct,
              var = ifelse(towct == 1, 0,
                           sum((n - mu)^2)/(towct - 1)),
              .groups = "drop") |>
    left_join(strata_wts, by = "strat") |>
    mutate(wt_mu = Area_SqNm * mu,
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
    group_by(sim, year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area,
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu,
              .groups = "drop") |>
    mutate(scenario = "Reallocation",
           pop = pop_num)
})

# Now group by pop (and sim if needed) to compute rel_ihat within each population
ihat_reall_all <- ihat_reall_all %>%
  group_by(pop, sim) %>%
  mutate(rel_ihat = stratmu / mean(stratmu)) %>%
  ungroup()


### Bind Indices ####
# bind all three scenario dataframes for efficient plotting and statistic calculation
#indices25 <- bind_rows(ihat_sq1, ihat_precl1, ihat_reall1)
indices25 <- bind_rows(ihat_sq_all,ihat_precl_all,ihat_reall_all)

# calculate the median relative abundance indices and the upper and lower confidence intervals across simulations and scenarios
# indices_med <- indices |>
#   group_by(year, scenario) |>
#   summarise(med = median(rel_ihat),
#             lower = quantile(rel_ihat, 0.025),
#             upper = quantile(rel_ihat, 0.975)) |>
#   mutate(type = "Estimated Relative Abundance Index")


## Plots ####

RelAbundPlot <- ggplot() +
  aes(x = year) +
  geom_line(data = trueN, aes(y = rel_N, color = scenario), linewidth = 1) +
  geom_point(data = indices25, aes(y = rel_ihat, color = fct_inorder(scenario)),
             position = position_dodge(width = 0.5)) +
  labs(x = "Year",
       y = "Relative Abundance (in kilograms)",
       subtitle = "Relative abundance of fall scup",
       color = "Scenario") +
  facet_wrap(~ pop, scales = "free_y", ncol = 1) +
  ylim(0, NA)



ggsave(str_c(species, season, "RelAbundPlot.png", sep = "_"),
       plot = RelAbundPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)

## SAVE THE DATA ####
saveRDS(trueN, here(surv.prods, str_c(species, season, "rel-TrueN.rds", sep = "_")))
saveRDS(indices25, here(surv.prods, str_c(species, season, "all-ihat-25survs.rds", sep = "_")))

saveRDS(ihat_sq_all, here(surv.prods, str_c(species, season, "25sims-sq_rel-ihat.rds", sep = "_")))
saveRDS(ihat_precl_all, here(surv.prods, str_c(species, season, "25sims-precl_rel-ihat.rds", sep = "_")))
saveRDS(ihat_reall_all, here(surv.prods, str_c(species, season, "25sims-reall_rel-ihat.rds", sep = "_")))
