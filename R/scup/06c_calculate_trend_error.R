### created: 01/22/2024
### updated: 05/01/2024

# 06c - CALCULATE ERROR IN TREND OVER TIME ####


## Objective ####
#

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(here)
library(broom)
# source(here("R", "sim_stratmean_fn.R"))


### DATA SET UP ####
#Directories
surv.prods <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")
plots <- here("outputs", "plots")

species <- "scup"
season <- "fall"
ages <- 0:7
years <- 1:15
nsims <- 1:100
nsurveys <- 25

# relative true abundance created here("R", "05_calculate_rel_abundance.R")
trueN <- readRDS(here(surv.prods, str_c(species, season, "all_rel-TrueN-100pops.rds", sep = "_")))

# relative abunance indices across scenarios created here("R", "05_calculate_rel_abundance.R")
indices <- readRDS(here(surv.prods, str_c(species, season, "all-ihat-25survs-100pops.rds", sep = "_")))


## ERROR IN TREND OVER TIME ####
# relative true abundance trends
trueN_lm <- trueN |>
  mutate(sim = 1) |>
  group_by(sim, scenario) |>
  nest() |>
  mutate(mods = map(data, ~lm(rel_N ~ year, data = .x)),
         coefs = map(mods, broom::tidy, conf.int = TRUE),
         slope = map(coefs, ~filter(.x, term == "year"))) |>
  select(sim, slope, scenario) |>
  unnest(cols = slope)

saveRDS(trueN_lm, here(surv.prods, str_c(species, season, "TrueRelAbundTrends.rds", sep = "_")))

indices_lm <- indices |>
  group_by(sim, scenario) |>
  nest() |>
  mutate(mods = map(data, ~lm(rel_ihat ~ year, data = .x)),
         coefs = map(mods, broom::tidy, conf.int = TRUE),
         slope = map(coefs, ~filter(.x, term == "year"))) |>
  select(sim, slope, scenario) |>
  unnest(cols = slope)

saveRDS(indices_lm, here(surv.prods, str_c(species, season, "RelIndexTrends.rds", sep = "_")))


slope_errors <- indices_lm |>
  select(sim, term, estimate) |>
  # left_join(trueN_lm[1:3], by = "sim") |>
  bind_cols(trueN_lm[1,1:3]) |>
  janitor::clean_names() |>
  mutate(err = estimate_4-estimate_7,
         rel_err = (err)/estimate_7,
         abs_rel_err = abs(rel_err))

saveRDS(slope_errors, here(surv.prods, str_c(species, season, "SlopeErrors.rds", sep = "_")))
#|>
  # group_by(scenario) |>
  # summarise(med_rel_err = median(rel_err),
  #           med_abs_err = median(abs_rel_err),
  #           lower_rel = quantile(rel_err, 0.025),
  #           upper_rel = quantile(rel_err, 0.975),
  #           lower_abs = quantile(abs_rel_err, 0.025),
  #           upper_abs = quantile(abs_rel_err, 0.975)) #|>
  #   ggplot() +
  # aes(x = as.factor(scenario), y = med_rel_err) +
  # geom_point(position = position_dodge2(width = 0.2)) +
  # geom_errorbar(aes(ymin=lower_rel, ymax=upper_rel), width= 0.2, position = position_dodge2(width = 0.2))


## PLOTS ####

# ggplot(indices_lm) +
#   geom_boxplot(aes(x = as.factor(scenario), y = estimate, color = scenario)) +
#   # ylim(0, NA) +
#   labs(x = "Year", y = "Linear regression slope estimates", title = str_c("Distribution of linear regression slope estimates for", season, species, "survey", sep = " ")) +
#   theme(legend.position = "bottom")

# Set the correct scenario order
indices_lm <- indices_lm |>
  mutate(scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation")))

slope_errors <- slope_errors |>
  mutate(scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation")))



#plots

ggplot(indices_lm) +
  geom_boxplot(aes(x = scenario, y = estimate, fill = fct_inorder(scenario)),
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") +
  scale_fill_manual(values = c("Status Quo" = "salmon", "Preclusion" = "goldenrod", "Reallocation" = "steelblue")) +
  labs(x = "Scenario", y = "Linear regression slope estimate",
       title = str_c("Distribution of linear regression slope estimates survey indices"), subtitle = "Scup - Fall") +
  theme_bw() +
  geom_hline(yintercept = trueN_lm$estimate[1], lty = 5, color = "chocolate4", linewidth = 0.75, show.legend = TRUE) +
  annotate("text", y = 0.01, x = 2, label = str_c("Relative true abundance slope: ", round(trueN_lm$estimate[1], 2)),
           color = "chocolate4", size = 5) +
  theme( text = element_text(size = 14),
         axis.title = element_text(size = 14),
         plot.title = element_text(size = 16, face = "bold"),
         legend.position = "none", legend.title = element_blank())



ggsave(str_c(species, season, "lm-ests-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)




ggplot(slope_errors) +
  geom_boxplot(aes(x = scenario, y = err, fill = fct_inorder(scenario)),
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") +
  scale_fill_manual(values = c("Status Quo" = "salmon", "Preclusion" = "goldenrod", "Reallocation" = "steelblue")) +
  labs(x = "Scenario", y = "Linear regression slope",
    title = str_c("Errors of linear regression slope estimates") , subtitle = "Scup - Fall") +
  theme_bw() +
  theme( text = element_text(size = 14),
         axis.title = element_text(size = 14),
         plot.title = element_text(size = 16, face = "bold"),
         legend.position = "none", legend.title = element_blank())


ggsave(str_c(species, season, "lm-errors-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)

ggplot(slope_errors) +
  geom_boxplot(aes(x = scenario, y = rel_err, fill = fct_inorder(scenario)),
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = "Relative errors of linear regression slope estimates",  subtitle = "Scup - Fall",
       x = "Scenario", y = "Linear regression slope", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               "Reallocation" = "steelblue")) +
  theme_bw() +
  theme( text = element_text(size = 14),
         axis.title = element_text(size = 14),
         plot.title = element_text(size = 16, face = "bold"),
         legend.position = "none", legend.title = element_blank())

ggsave(str_c(species, season, "lm-rel-errors-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)


ggplot(slope_errors) +
  geom_boxplot(aes(x = scenario, y = abs_rel_err, fill = fct_inorder(scenario)),
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") +
  labs(title = "Absolute relative errors of linear regression slope estimates", subtitle = "Scup - Fall",
       x = "Scenario", y = "Linear regression slope", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               "Reallocation" = "steelblue")) +
  theme_bw() +
  theme( text = element_text(size = 14),
         axis.title = element_text(size = 14),
         plot.title = element_text(size = 16, face = "bold"),
         legend.position = "none", legend.title = element_blank())

ggsave(str_c(species, season, "lm-abs-rel-errors-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 9, height = 6)

