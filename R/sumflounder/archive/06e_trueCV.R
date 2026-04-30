
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
species <- "summerflounder"
season  <- "spring"
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
  ylim(0, 1.5) +
  labs(title = " ", subtitle = "Scup - Fall",
       x = "Year", y = "True CV", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo"   = "#D76E9A",
                               "Preclusion"   = "#E9988C",
                               "Reallocation" = "#785838",
                               "Ratio estimator" = "#426737",
                               "Model based" = "#93995C",
                               "Model based wind" = "#FDE16A")) +
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
    position = position_dodge(width = 0.8),     outlier.shape = NA, color = "grey30") +
  geom_line(data = true_N_cv, aes(x = factor(year), y = cv, group = 1), color = "grey10",linewidth = 1) +
  geom_point(data = true_N_cv, aes(x = factor(year), y = cv), color = "grey20",size = 3) +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(title = "", subtitle = "Scup - Fall", x = "Year", y = "CV", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo"   = "#D76E9A",
                               "Preclusion"   = "#E9988C",
                               "Reallocation" = "#785838",
                               "Ratio estimator" = "#426737",
                               "Model based" = "#93995C",
                               "Model based wind" = "#FDE16A")) +
  theme_bw() + theme(text = element_text(size = 20), axis.title = element_text(size = 20),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right", legend.title = element_blank())

TRUECVboxplot_merged

ggsave(str_c(species, season, "TRUECVboxplot_merged.png", sep = "_"),
       plot = TRUECVboxplot_merged,
       device = "png",
       # last_plot(),
       here(plots),
       width = 10, height = 6)
