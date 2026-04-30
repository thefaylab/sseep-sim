
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
indices <- readRDS(here(surv.prods, str_c(species, season, "indices.rds", sep="_")))



true_cv <- indices %>%
  group_by(pop, year, scenario) %>%
  summarise(
    true_cv = sd(rel_ihat, na.rm = TRUE) / mean(rel_ihat, na.rm = TRUE),
    .groups = "drop"
  )

true_cv_plot <- true_cv %>%
  mutate(
    scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation",
                                           "Ratio estimator", "Model based", "Model based wind"))
  )

TRUECVboxplot <- ggplot(true_cv_plot, aes(x= factor(year), y = true_cv, fill = scenario)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  ylim(0, 0.3) +
  labs(title = " ", subtitle = "Scup - Fall",
       x = "Year", y = "True CV", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo"   = "#2F397A",
                               "Preclusion"   = "#7391BD",
                               "Reallocation" = "#894846",
                               "Ratio estimator" = "#785838",
                               "Model based" = "#93995C",
                               "Model based wind" = "#4F6009")) +
  theme(text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "right", legend.title = element_blank())



true_N_cv <- trueN %>%
  group_by(year) %>%
  summarise(
    cv = sd(rel_N, na.rm = TRUE) / mean(rel_N, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(true_N_cv, aes(x= factor(year), y = cv)) +
  geom_point(color = "black") +
  ylim(0, 0.7) +
  labs(title = "", subtitle = "Scup - Fall",
    x = "Year", y = "CV true rel_N") +
  theme_bw() +
  theme(text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )


TRUECVboxplot_merged <- ggplot() +
  geom_boxplot(data = true_cv_plot, aes(x = factor(year), y = true_cv, fill = scenario),
    position = position_dodge(width = 0.8),     outlier.shape = NA, color = "black") +
  geom_line(data = true_N_cv, aes(x = factor(year), y = cv, group = 1), color = "black",linewidth = 1.2) +
  geom_point(data = true_N_cv, aes(x = factor(year), y = cv), color = "black",size = 3) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(title = "", subtitle = "Scup - Fall", x = "Year", y = "CV", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "#2F397A",
                               "Preclusion" = "#7391BD",
                               "Reallocation" = "#894846",
                               "Ratio estimator" = "#785838",
                               "Model based" = "#93995C",
                               "Model based wind" = "#4F6009")) +
  theme_bw() + theme(text = element_text(size = 20), axis.title = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right", legend.title = element_blank())

TRUECVboxplot_merged
