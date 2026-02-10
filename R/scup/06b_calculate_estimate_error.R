### created: 01/18/2024
### updated: 02/07/2025

# 06b - CALCULATE ESTIMATION ERROR ####


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
errors <- readRDS(here(perform.metrics, str_c(species, season,"25all-rel-error-100pops.rds", sep = "_")))


## ESTIMATION ERROR ####
errors <- errors %>%
  mutate(scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation")),
         period = if_else(year <= 5, "Years 1-5", "Years 6-15"))


CVboxplot <- ggplot(errors, aes(x = period, y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") + ylim(0,1) +
  labs(title = "Distribution of CVs ",
       x = "Period", y = "Relative Error", fill = "Scenario") +
  scale_fill_manual(values = c("Status Quo" = "salmon",
                               "Preclusion" = "goldenrod",
                               "Reallocation" = "steelblue")) +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right", legend.title = element_blank())




ggsave(str_c(species, season, "CV-boxplot.png", sep = "_"),
       plot = CVboxplot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)
















ggplot(errors, aes(x = scenario, y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = 21, outlier.fill = "white", color = "black") +
  labs(
    title = "Distribution of Coefficients of Variation",
    x = "Scenario",
    y = "CV",
    fill = "Scenario"
  ) +
  scale_fill_manual(values = c("Status Quo" = "salmon",  # Default R Orange
                               "Preclusion" = "#DAA520",  # Default R Green
                               "Reallocation" = "#4682B4")) + # Default R Blue
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none")



# First filter to only years 6 to 15
errors_filtered <- errors %>%
  filter(year >= 6)

# Then plot by individual years (6–15)
ggplot(errors_filtered, aes(x = as.factor(year), y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.shape = NA,
               color = "black") +
  ylim(0, .75) +
  labs(title = "Distribution of CVs (Years 6–15)",
       x = "Year", y = "Coefficient of Variation", fill = "Scenario") +
  scale_fill_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue"
  )) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.title = element_blank()
  )




pops_to_plot <- c(2,5,88)
sim_to_plot <- c(1,15)


library(tidyverse)
indices_cv_wide <- indices %>%
  select(sim, year, pop, scenario, cv) %>%  # select only needed columns
  pivot_wider(names_from = scenario, values_from = cv)


indices_cv_wide <- indices_cv_wide %>%
  mutate(
    diff_SQ_Preclusion = `Status Quo` - Preclusion,
    diff_SQ_Reallocation = `Status Quo` - Reallocation
  )


indices_cv_wide_sub <- filter(indices_cv_wide, pop %in% pops_to_plot, sim %in% sim_to_plot)


indices_cv_wide_sub <- indices_cv_wide_sub %>%
  left_join(trueN_yr %>% select(year, cv_N), by = "year")

indices_cv_wide_sub <- indices_cv_wide_sub %>%
  mutate(
    diff_true_SQ = `Status Quo` - cv_N,
    diff_true_Preclusion = Preclusion - cv_N,
    diff_true_Reallocation = Reallocation - cv_N
  )

custom_colors <- c(
  "Status Quo" = "salmon",
  "Preclusion" = "goldenrod",
  "Reallocation" = "steelblue",
  "cv_N" = "grey60"
)


indices_cv_long_sub <- indices_cv_wide_sub %>%
  pivot_longer(cols = c(`Status Quo`, Preclusion, Reallocation, cv_N),
               names_to = "Scenario", values_to = "cv")

indices_cv_long_sub <- indices_cv_long_sub %>%
  mutate(
    Scenario = factor(Scenario, levels = c("cv_N","Status Quo", "Preclusion", "Reallocation"))
  )


ggplot(indices_cv_long_sub, aes(x = year, y = cv, color = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(sim ~ pop) +   # switch strips to top and right
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Comparison of Survey and True CVs by Population and Simulation",
    x = "Year",
    y = "Coefficient of Variation (CV)",
    color = "Scenario"
  ) +
  theme_minimal() +
  theme(
    strip.text.x.top = element_text(size = 10, face = "bold"),   # top facet (pop)
    strip.text.y.right = element_text(size = 10, face = "bold"), # right facet (sim)
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    plot.margin = margin(10, 30, 10, 30)
  )




# Short and simple period labels
indices_cv_wide_sub <- indices_cv_wide_sub %>%
  mutate(
    period = case_when(
      year %in% 1:5 ~ "Years 1–5",
      year %in% 6:15 ~ "Years 6–15"
    )
  )

# Pivot longer
indices_cv_long_sub <- indices_cv_wide_sub %>%
  pivot_longer(cols = c(`Status Quo`, Preclusion, Reallocation),
               names_to = "scenario", values_to = "cv")





# Make sure 'scenario' has correct order
indices_cv_long_sub <- indices_cv_long_sub %>%
  mutate(
    scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation"))
  )

# Create new LABELLER:
custom_labeller <- labeller(
  pop = c("2" = "Population 2", "5" = "Population 5", "88" = "Population 88"),
  sim = c("1" = "Survey 1", "15" ="Survey 15")
)

# Plot
ggplot(indices_cv_long_sub, aes(x = period, y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  facet_grid(sim ~ pop, labeller = custom_labeller) +
  # <-- key: use sim ~ pop to preserve your layout
  scale_fill_manual(
    values = c(
      "Status Quo" = "salmon",
      "Preclusion" = "goldenrod",
      "Reallocation" = "steelblue"
    ),
    breaks = c("Status Quo", "Preclusion", "Reallocation")
  ) +
  labs(
    title = "Subset of survey CVs by population and year periods",
    x = "Period",
    y = "Coefficient of Variation (CV)",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    strip.text.x.top = element_text(size = 10, face = "bold"),    # Top facet text (Population)
    strip.text.y.right = element_text(size = 10, face = "bold"),  # Right facet text (Survey)
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    plot.margin = margin(10, 30, 10, 30)
  )


# Add a new "scenario" for True N CV
trueN_rows <- indices_cv_long_sub %>%
  mutate(
    scenario = "True N CV",
    cv = cv_N  # set cv equal to cv_N
  )

# Combine with original
indices_cv_long_sub2 <- bind_rows(
  indices_cv_long_sub,
  trueN_rows
)


ggplot(indices_cv_long_sub2, aes(x = period, y = cv, fill = scenario)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  facet_grid(sim ~ pop, labeller = custom_labeller) +
  scale_fill_manual(
    values = c(
      "Status Quo" = "salmon",
      "Preclusion" = "goldenrod",
      "Reallocation" = "steelblue",
      "True N CV" = "grey"
    ),
    breaks = c("Status Quo", "Preclusion", "Reallocation", "True N CV")
  ) +
  labs(
    title = "Subset of survey and true CVs by population and year periods",
    x = "Period",
    y = "Coefficient of Variation (CV)",
    fill = "Scenario"
  ) +
  theme_minimal() +
  theme(
    strip.text.x.top = element_text(size = 10, face = "bold"),
    strip.text.y.right = element_text(size = 10, face = "bold"),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing = unit(0.3, "lines"),
    plot.margin = margin(10, 30, 10, 30)
  )


