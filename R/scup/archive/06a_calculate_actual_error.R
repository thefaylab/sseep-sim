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


#Indices

# relative true abundance created here("R", "05_calculate_rel_abundance.R")
trueN <- readRDS(here(surv.prods, str_c(species, season, "TrueN.rds", sep = "_")))

# relative abunance indices across scenarios created here("R", "05_calculate_rel_abundance.R")
indices <- readRDS(here(surv.prods, str_c(species, season, "indices.rds", sep = "_")))


## CALCULATE RELATIVE AND ABSOLUTE ERRORS ####
errors <- indices %>%
  mutate(sim = as.character(sim), pop = as.character(pop)) |>
  left_join( trueN |> mutate(pop = as.character(pop)), by = c("pop", "year")) |>
  mutate(rel_err = (rel_ihat - rel_N) / rel_N,
         abs_rel_err = abs(rel_err)) |>
  select(!scenario.y) |> #delete TRUE name column only.
  rename(scenario = scenario.x)#scenarios and populations are relative to TRUE rel N


## PLOTS ####
# relative error plot

errors_plot <- errors %>%
  filter((scenario == "Status Quo") | (scenario %in% c("Preclusion", "Reallocation",
                                                       "Ratio estimator", "Model based","Model based wind") & year >= 6)) |>
  mutate(scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation",
                                                "Ratio estimator", "Model based", "Model based wind")))



RelErrBoxPlot <- ggplot(errors_plot, aes(x = , y = rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(-1.2,1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#535260", linewidth = 1) +
  labs(title = "", subtitle = "Scup - Fall",
       x = NULL, y = "Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo"   = "#2F397A",
                               "Preclusion"   = "#7391BD",
                               "Reallocation" = "#894846",
                               "Ratio estimator" = "#785838",
                               "Model based" = "#93995C",
                               "Model based wind" = "#4F6009")) +
  theme(text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "right", legend.title = element_blank())


ggsave(str_c(species, season, "RelErrBoxPlot.png", sep = "_"),
       plot = RelErrBoxPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)



##Absolute Error

AbsRelErrBoxPlot <- ggplot(errors_plot, aes(x = , y = rel_err, fill = scenario))  +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(0,3) +
  #  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = "Distribution of Absolute Relative Errors", subtitle = "Scup - Fall",
       x = NULL, y = "Absolute Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo"   = "#2F397A",
                               "Preclusion"   = "#7391BD",
                               "Reallocation" = "#894846",
                               "Ratio estimator" = "#785838",
                               "Model based" = "#93995C",
                               "Model based wind" = "#4F6009")) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right", legend.title = element_blank())


ggsave(str_c(species, season, "AbsRelErrBoxPlot.png", sep = "_"),
       plot = AbsRelErrBoxPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)




## SAVE THE DATA ####
saveRDS(errors, here(perform.metrics, str_c(species, season, "rel-error.rds", sep = "_")))



