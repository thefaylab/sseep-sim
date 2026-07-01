
### PACKAGES ####
library(tidyverse)
library(here)
library(sdmTMB)
library(SimSurvey)
source(here("R", "sim_stratmean_fn.R"))

# estimates must have:
# pop, sim, year, scenario, ihat, cv
# where:
# ihat = estimated abundance index or relative index
# cv   = perceived/estimated CV from the variance formula or model SE

### DATA SET UP ####
# Directories
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

# Parameters
species <- "scup"
season  <- "fall"
ages      <- 0:7
years     <- 1:15
nsims   <- 1:100
ids     <- sprintf("%03d", nsims)


# Load true population and estimated indices
trueN <- readRDS(here(surv.prods, str_c(species, season, "trueN.rds", sep="_")))

indices <- readRDS(here(surv.prods, str_c(species, season, "indices.rds", sep="_"))) |>
  select(sim, year, stratmu, stratvar,  cv, scenario, pop, est)


cmp <- indices |>
  mutate(index_value = if_else(
    scenario %in% c("Model based", "Model based wind"),
    est,
    stratmu)) |>
  group_by(pop, year, scenario) |>
  mutate(true_cv = sd(index_value, na.rm = TRUE) / mean(index_value, na.rm = TRUE)) |>
  ungroup()


# cmp_long <- cmp |>
#   select(pop, sim, scenario, year, cv, true_cv) |>
#   pivot_longer(cols = c(cv, true_cv), names_to = "cv_type", values_to = "cv_value") |>
#   mutate(cv_type = recode(cv_type, cv = "Estimated CV", true_cv = "True CV"),
#          fill_group = if_else(cv_type == "True CV", "True CV", scenario)) |>
#     mutate(scenario = factor(scenario,
#     levels = c("Status Quo", "Preclusion", "Reallocation",
#                "Ratio estimator", "Model based", "Model based wind")))
#
# true_cv <- ggplot(cmp_long, aes(x = factor(year), y = cv_value, fill = fill_group)) +
#   geom_boxplot(aes(group = interaction(year, cv_type)),
#                position = position_dodge(0.75), width = 0.6, outlier.shape = NA) +
#   facet_wrap(~ scenario, scales = "free_y") +
#   scale_fill_manual(values = c(
#     "Status Quo" = "#9B5F6B",
#     "Preclusion" = "#8C2B0E",
#     "Reallocation" = "#C5692D",
#     "Ratio estimator" = "#D9A441",
#     "Model based" = "#A7C957",
#     "Model based wind" = "#275E4D",
#     "True CV" = "#D4D9DD")) +
#   labs(x = "Year", y = "CV", title = "Estimated vs True CV",
#        subtitle = "Scup - Fall",fill = NULL) +
#   theme_bw() +
# theme(text = element_text(size = 14),
#       axis.title = element_text(size = 14),
#       plot.title = element_text(size = 16, face = "bold"),
#       legend.position = "right", legend.title = element_blank())
#
#
# ggsave(str_c(species, season, "true_cv.png", sep = "_"),
#        plot = true_cv,
#        device = "png",
#        # last_plot(),
#        here(plots),
#        width = 15, height = 10)
#
#
# cmp_diff <- cmp |>
#   select(pop, sim, scenario, year, cv, true_cv) |>
#   mutate(cv_diff = cv - true_cv,
#          scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation",
#                                                 "Ratio estimator", "Model based", "Model based wind"))) |>
#   arrange(scenario, pop, sim, year)
#
#
#
# cv_diff <- ggplot(cmp_diff, aes(x = factor(year), y = cv_diff, fill = scenario)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_boxplot(outlier.shape = NA, width = 0.6) +
#   facet_wrap(~ scenario, scales = "free_y") +
#   scale_fill_manual(values = c("Status Quo" = "#9B5F6B",
#                                "Preclusion" = "#8C2B0E",
#                                "Reallocation" = "#C5692D",
#                                "Ratio estimator" = "#D9A441",
#                                "Model based" = "#A7C957",
#                                "Model based wind" = "#275E4D")) +
#   labs(x = "Year", y = "True CV - Estimated CV",
#        title = "Difference between True and Estimated CV",
#        subtitle = "Scup - Fall", fill = NULL) +
#   theme_bw() +
#   theme(text = element_text(size = 14), axis.title = element_text(size = 14),
#         plot.title = element_text(size = 16, face = "bold"),
#         legend.position = "none")
#
# ggsave(str_c(species, season, "cv_diff.png", sep = "_"),
#        plot = cv_diff,
#        device = "png",
#        # last_plot(),
#        here(plots),
#        width = 15, height = 10)


cmp_diff <- cmp |>
  select(pop, sim, scenario, year, cv, true_cv) |>
  filter(year >= 6) |>
  mutate(cv_diff = cv - true_cv,
         scenario = factor(scenario,
                           levels = c("Status Quo", "Preclusion",
                                      "Reallocation", "Ratio estimator",
                                      "Model based", "Model based wind"))) |>
  arrange(scenario, pop, sim, year)

cv_diff <- ggplot(cmp_diff, aes(x = , y = cv_diff, fill = scenario))  +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  ylim(-3,1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "chocolate4", linewidth = 1) +
  labs(title = "Difference between estimated and true CV", subtitle = "Scup - Fall",
       x = NULL, y = "Est - True cv", fill = "NULL") +
  scale_fill_manual(values = c(
    "Status Quo"       = "#9B5F6B",
    "Preclusion"       = "#8C2B0E",
    "Reallocation"     = "#C5692D",
    "Ratio estimator"  = "#D9A441",
    "Model based"      = "#A7C957",
    "Model based wind" = "#275E4D")) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right", legend.title = element_blank())

ggsave(str_c(species, season, "cv_diff.png", sep = "_"),
       plot = cv_diff,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)


