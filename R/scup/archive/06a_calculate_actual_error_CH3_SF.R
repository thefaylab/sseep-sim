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
library(patchwork)
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
# survdat_sq <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, .x))))
# survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))
# survdat_reall <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_reall_survey.rds", species, season, .x))))
# dist          <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))



# relative true abundance created here("R", "05_calculate_rel_abundance.R")
trueN <- readRDS(here(surv.prods, str_c(species, season, "rel-TrueN-100pops.rds", sep = "_")))

# relative abunance indices across scenarios created here("R", "05_calculate_rel_abundance.R")
indices_fixed <- readRDS(here(surv.prods, str_c(species, season, "all-ihat-25survs-100pops-hybrids.rds", sep = "_")))

#i1to5_pop1 <- indices |> filter(year <= 5 & pop == 1 & year == 1 & sim <= 5)

## CALCULATE RELATIVE AND ABSOLUTE ERRORS ####
errors_fixed <- indices_fixed %>%
  mutate(sim = as.character(sim), pop = as.character(pop)) |>
  left_join( trueN |> mutate(pop = as.character(pop)), by = c("pop", "year")) |>
  mutate(rel_err = (rel_ihat - rel_N) / rel_N,
         abs_rel_err = abs(rel_err)) |>
  select(!scenario.y) |> #delete TRUE name column only.
  rename(scenario = scenario.x)  #scenarios and populations are relative to TRUE rel N



## PLOTS ####
# relative error plot

errors_plot_s2 <- errors_fixed %>%
  mutate(
    period = case_when(
      year <= 5 ~ "Years 1–5",
      year >= 6 ~ "Years 6–15"
    )
  ) %>%
  filter((scenario %in% c("Status Quo","Preclusion", "Hybrid") & year >= 6)) %>%
  mutate(
    scenario = dplyr::recode(scenario, "Hybrid" = "Supplemental"),
    scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Supplemental"))
  )



RelErrBoxPlot_year_s2 <- ggplot(errors_plot_s2, aes(x = as.factor(year), y = rel_err)) +
  geom_boxplot(aes(fill = scenario),
               position = position_dodge2(preserve = "single", width = 0.8),
               width = 0.7,
               outlier.shape = NA,
               color = "black") + ylim(-3,3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                              # "Hybrid Small" = "seagreen",
                               "Hybrid" = "slateblue")) +
  labs(x = "Year", y = "Relative Error", fill = "Scenario",
       subtitle = "Scup - Fall") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.title = element_blank()
  )





RelErrBoxPlot_period_s2 <- ggplot(errors_plot_s2, aes(y = rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge2(preserve = "single", width = 0.4),
               width = 0.4,
               outlier.shape = NA,
               color = "black") + ylim(-1,1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = " ", subtitle = "Scup",
       x = "", y = "Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               # "Hybrid Small" = "seagreen",
                               "Supplemental" = "slateblue")) +
  theme_bw() +
  theme(text = element_text(size = 30),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ggsave(str_c(species, season, "RelErrBoxPlot_period_s2.png", sep = "_"),
       plot = RelErrBoxPlot_period_s2,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 7)



##Absolute Error

AbsRelErrBoxPlot <- ggplot(errors_plot_s2, aes(y = abs_rel_err, fill = scenario)) +
  geom_boxplot(position = position_dodge2(preserve = "single", width = 0.9),
               width = 0.7,
               outlier.shape = NA,
               color = "black") + ylim(0,1.15) +
  labs(title = " ", subtitle = "Scup",
       x = "", y = "Absolute Rel Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               # "Hybrid Small" = "seagreen",
                               "Supplemental" = "slateblue")) +
  theme_bw() +
  theme(text = element_text(size = 30),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(str_c(species, season, "AbsRelErrBoxPlot.png", sep = "_"),
       plot = AbsRelErrBoxPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 7)



CVboxplot<- ggplot(errors_plot_s2, aes(y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge2(preserve = "single", width = 0.9),
               width = 0.7,
               outlier.shape = NA,
               color = "black") + ylim(0,0.6) +
  labs(title = " ", subtitle = "Scup",
       x = "", y = "CV", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                              # "Hybrid Small" = "seagreen",
                               "Supplemental" = "slateblue")) +
  theme_bw() +
  theme(text = element_text(size = 30),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(str_c(species, season, "CVboxplot.png", sep = "_"),
       plot = CVboxplot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 7)



RelErrBoxPlot_period_s2 + AbsRelErrBoxPlot + CVboxplot






errors_k <- ihat_hybrid_sens_all %>%
  mutate(sim = as.character(sim), pop = as.character(pop)) |>
  left_join( trueN |> mutate(pop = as.character(pop)), by = c("pop", "year")) |>
  mutate(rel_err = (rel_ihat - rel_N) / rel_N,
         abs_rel_err = abs(rel_err))  |>
  select(!scenario.y) |> #delete TRUE name column only.
  rename(scenario = scenario.x)#scenarios and populations are relative to TRUE rel N


errors_k2 <- errors_k %>%
  mutate(
    period = case_when(
      year <= 5 ~ "Years 1–5",
      year >= 6 ~ "Years 6–15"
    )
  )


ggplot(errors_k2 |> filter(year >= 6), aes(x = period, y = rel_err)) +
  geom_boxplot(aes(fill = factor(k)),
               position = position_dodge2(preserve = "single", width = 0.9),
               width = 0.7,
               outlier.shape = NA,
               color = "black") + ylim(-1.5,1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(x = "", y = "Relative Error", fill = "k",
       subtitle = "Scup") +
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold"),
        legend.position = "right", legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




ggplot(errors_k2 |> filter(year >= 6), aes(y = cv, fill = factor(k))) +
  geom_boxplot(position = position_dodge2(preserve = "single", width = 0.4),
               width = 0.7,
               outlier.shape = NA,
               color = "black") + ylim(0,0.2) +
  labs(title = " ", subtitle = "Scup",
       x = "", y = "CV", fill = "k") +
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold"),
        legend.position = "right", legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
