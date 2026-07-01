### created: 01/18/2024
### updated: 02/07/2025

# 06b - CALCULATE PERCEIVED UNCERTAINTY ####


## Objective ####
# For a given species, distribution, and survey, calulate


### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
# source(here("R", "sim_stratmean_fn.R"))


### DATA SET UP ####
#Directories
surv.prod <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")
plots <- here("outputs", "plots")

species <- "scup"
season <- "fall"
ages <- 0:7
years <- 1:15
nsims <- 1:100
nsurveys <- 25

### LOAD DATA ####
# relative and absolute errors of simulated abundance created here("R", "06a_calculate_actual_abundance.R")
errors <- readRDS(here(perform.metrics, str_c(species, season,"rel-error.rds", sep = "_")))

errors |> filter(sim==1, year==1, pop==1)
## ESTIMATION ERROR ####
errors_plot <- errors %>%
  filter((scenario == "Status Quo") |
           (scenario %in% c("Preclusion", "Supplemental Survey", "Supp correction") & year >= 6)) |>
  mutate(scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Supplemental Survey", "Supp correction")))



CVboxplot_Ch3<- ggplot(errors_plot, aes(y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(0,1) +
  labs(title = " ", subtitle = "Scup - Fall",
       x = "", y = "CV", fill = "Scenario") +
  scale_fill_manual(values = c(
    "Status Quo"       = "#9B5F6B",
    "Preclusion"       = "#8C2B0E",
    "Supplemental Survey"     = "#6592B0",
    "Supp correction"  = "#0067A2"
  )) +
  theme(text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "right", legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



ggsave(str_c(species, season, "CV-boxplot_Ch3.png", sep = "_"),
       plot = CVboxplot_Ch3,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)


