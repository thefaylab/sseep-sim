### created: 02/05/2024
### updated: 02/07/2025

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
#Directories
surv.prods <- here("data", "rds", "surv-prods")
perform.metrics <- here("data", "rds", "perform-metrics")
plots <- here("outputs", "plots")

species <- "summerflounder"
season <- "spring"
ages <- 0:7
years <- 1:15
nsims <- 1:100
nsurveys <- 25

### LOAD DATA ####
# relative abunance indices for each survey effort scenario created here("R", "05_calculate_rel_abundance.R")
### STATUS QUO
ihat_sq_all <- readRDS(here(surv.prods, str_c(species, season, "100pops-25sims-sq_rel-ihat.rds", sep = "_")))

### PRECLUSION
ihat_precl_all <- readRDS(here(surv.prods, str_c(species, season, "100pops-25sims-precl_rel-ihat.rds", sep = "_")))

### REALLOCATION
ihat_reall_all <- readRDS(here(surv.prods, str_c(species, season, "100pops-25sims-reall_rel-ihat.rds", sep = "_")))

### errors
errors <- readRDS(here(perform.metrics, str_c(species, season,"25all-rel-error-100pops.rds", sep = "_")))


## CALCULATE ABSOLUTE PERCENT DIFFERENCES ####
### Status Quo vs Preclusion ####
sq_precl_diff <- left_join(ihat_sq_all, ihat_precl_all, by = c("pop","sim", "year")) |>
  janitor::clean_names() |>
  group_by(pop,sim, year) |>
  mutate(perc_diff_mu =( abs(rel_ihat_x - rel_ihat_y) / rel_ihat_x)*100,
         perc_diff_cv =( abs(cv_x - cv_y) / cv_x)*100) |>
  select(pop,sim, year, perc_diff_mu, perc_diff_cv) |>
  mutate(type = "SQ-Preclusion")



### Status Quo vs Reallocation ####
sq_reall_diff <- left_join(ihat_sq_all, ihat_reall_all, by = c("pop","sim", "year")) |>
  janitor::clean_names() |>
  group_by(pop,sim, year) |>
  mutate(perc_diff_mu =( abs(rel_ihat_x - rel_ihat_y) / rel_ihat_x)*100,
         perc_diff_cv =( abs(cv_x - cv_y) / cv_x)*100) |>
  select(pop,sim, year, perc_diff_mu, perc_diff_cv) |>
  mutate(type = "SQ-Reallocation")


## BIND ####
diffs <- bind_rows(sq_precl_diff, sq_reall_diff)




## CALCULATE ABSOLUTE PERCENT DIFFERENCES for relative errors ####
### Status Quo vs Preclusion ####
sq_precl_relerr <- left_join(errors %>% filter(scenario == "Status Quo"),
                             errors %>% filter(scenario == "Preclusion"),
                             by = c("pop", "sim", "year")) |>
  janitor::clean_names() |>
  group_by(pop, sim, year) |>
  mutate(diff_relerr = (rel_err_x - rel_err_y),
         diff_abserr = abs(rel_err_x - rel_err_y)) |>
  select(pop, sim, year, diff_relerr, diff_abserr) |>
  mutate(type = "SQ-Preclusion")

### Status Quo vs Reallocation ####
sq_reall_relerr <- left_join(errors %>% filter(scenario == "Status Quo"),
                             errors %>% filter(scenario == "Reallocation"),
                             by = c("pop", "sim", "year")) |>
  janitor::clean_names() |>
  group_by(pop, sim, year) |>
  mutate(diff_relerr = (rel_err_x - rel_err_y),
         diff_abserr = abs(rel_err_x - rel_err_y)) |>
  select(pop, sim, year, diff_relerr, diff_abserr) |>
  mutate(type = "SQ-Reallocation")

## BIND ####
relerr_diffs <- bind_rows(sq_precl_relerr, sq_reall_relerr)

relerr_diffs <- relerr_diffs %>%
  mutate(period = ifelse(year <= 5, "1-5", "6-15"))
relerr_diffs_late <- relerr_diffs %>%
  filter(period == "6-15")

ggplot(relerr_diffs, aes(x = type, y = diff_relerr, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.75),
               outlier.shape = NA, color = "black") +
  facet_wrap(~ period) +
  scale_fill_manual(values = c("SQ-Preclusion" = "goldenrod4",
                               "SQ-Reallocation" = "steelblue4")) +
  ylim(0,1) +
  labs(title = "Difference in relative error among scenarios and time period",
       subtitle = "Summerflounder - Spring",
       x = "Scenario",
       y = "Relative Error Difference",
       fill = "Scenario") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14),
        legend.position = "none")


ggsave(str_c(species, season, "relerr-diffs-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)



ggplot(relerr_diffs_late, aes(x = factor(year), y = diff_relerr, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.75),
               outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("SQ-Preclusion" = "goldenrod4",
                               "SQ-Reallocation" = "steelblue4")) +
  ylim(-1.5, 1.5) +
  labs(title = "Difference in relative error among scenarios",
       subtitle = "Summerflounder - Spring",
       x = "Year",
       y = "Relative error difference",
       fill = "Scenario") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom")

ggsave(str_c(species, season, "relerr-diffs-late-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)



"Difference in absolute relative error among scenarios (Years 6-15)"
ggplot(relerr_diffs_late, aes(x = factor(year), y = diff_abserr, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.75),
               outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("SQ-Preclusion" = "goldenrod4",
                               "SQ-Reallocation" = "steelblue4")) +
  ylim(0, 1) +
  labs(title = "Difference in absolute relative error among scenario",
       subtitle = "Summerflounder - Spring",
       x = "Year",
       y = "Absolute relative error difference",
       fill = "Scenario") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom")

ggsave(str_c(species, season, "absrelerr-diffs-late-boxplot.png", sep = "_"), device = "png", last_plot(), here(plots), width = 8, height = 6)
