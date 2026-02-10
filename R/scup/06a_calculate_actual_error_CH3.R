### created: 01/18/2024
### updated: 02/07/2025

# 06a - CALCULATE ACTUAL ERROR ####


## Objective ####
# For a given species and iteration, calculate the relative and absolute error to compare an abundance index to the true abundance

# Outputs:
#  a dataframe with values for each simulation and year pertaining to the relative error and absolute relative error of a given survey and its abundance index compared to the relative true abundance
# the relative errors and absolute errors for a given survey scenario plotted across time


### PACKAGES ####
library(tidyverse)
library(here)
library(data.table)
library(sdmTMB)
source(here("R", "sim_stratmean_fn.R"))
theme_set(theme_bw())


### DATA SET UP ####
#Directories
surv.prods <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")
plots <- here("outputs", "plots")

# Parameters
species <- "scup"
season  <- "fall"
ages      <- 0:7
years     <- 1:15
nsims   <- 1:100
ids     <- sprintf("%03d", nsims)
nsurveys <- 25


#Data
# pop <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_abund-dist.rds", species, season, .x))))
 survdat_sq <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, .x))))
# survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))
# survdat_reall <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_reall_survey.rds", species, season, .x))))
# dist          <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))



# relative true abundance created here("R", "05_calculate_rel_abundance.R")
trueN <- readRDS(here(surv.prods, str_c(species, season, "all_rel-TrueN-100pops.rds", sep = "_")))

# relative abunance indices across scenarios created here("R", "05_calculate_rel_abundance.R")
indices_2 <- readRDS(here(surv.prods, str_c(species, season, "all-ihat-25survs-100pops.rds", sep = "_")))

#i1to5_pop1 <- indices |> filter(year <= 5 & pop == 1 & year == 1 & sim <= 5)

## CALCULATE RELATIVE AND ABSOLUTE ERRORS ####
errors <- indices2 %>%
  mutate(sim = as.character(sim), pop = as.character(pop)) |>
  left_join( trueN |> mutate(pop = as.character(pop)), by = c("pop", "year")) |>
  mutate(rel_err = (rel_ihat - rel_N) / rel_N,
         abs_rel_err = abs(rel_err)) |>
  select(!scenario.y) |> #delete TRUE name column only.
  rename(scenario = scenario.x)#scenarios and populations are relative to TRUE rel N


## PLOTS ####
# relative error plot

errors_plot <- errors %>%
  filter(
    (scenario == "Status Quo") |
      (scenario %in% c("Preclusion", "Hybrid") & year >= 6)
  )


RelErrBoxPlot <- ggplot(errors, aes(x = , y = rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(-1,1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = "", subtitle = "Scup - Fall",
       x = "Year", y = "Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               "Hybrid" = "seagreen")) +
  theme(text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "right", legend.title = element_blank())


ggsave(str_c(species, season, "RelErrBoxPlot.png", sep = "_"),
       plot = RelErrBoxPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 7, height = 5)



##Absolute Error

AbsRelErrBoxPlot <- ggplot(errors, aes(x = period, y = abs_rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(0,1.5) +
  #  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = "Distribution of Absolute Relative Errors", subtitle = "Scup - Fall",
       x = "Period", y = "Absolute Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               "Reallocation" = "steelblue",
                               "Hybrid" = "seagreen")) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right", legend.title = element_blank())


ggsave(str_c(species, season, "AbsRelErrBoxPlot.png", sep = "_"),
       plot = AbsRelErrBoxPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)




## SAVE THE DATA ####
saveRDS(errors, here(perform.metrics, str_c(species, season, "25all-rel-error-100pops.rds", sep = "_")))

























##RElative Error

ggplot(errors, aes(x = as.factor(year), y = rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(-1.5,1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = "Distribution of Relative Errors",
    x = "Year", y = "Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               "Reallocation" = "steelblue")) +
  theme( text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right", legend.title = element_blank())






ggplot(errors, aes(x = scenario, y = rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = 21, outlier.fill = "white", color = "black") +
  labs(
    title = "Distribution of Relative Errors",
    x = "Scenario",
    y = "Relative Error",
    fill = "Scenario"
  ) +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "#DAA520",
                               "Reallocation" = "#4682B4")) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none")






