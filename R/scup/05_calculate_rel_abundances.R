### created: 01/18/2024
### updated: 02/07/2024

# 05 - CALCULATE RELATIVE ABUNDANCES ####


## Objective ####
# For a given species, distribution, and survey, calculate the relative true abundance and the relative abundance index.

# Outputs: Relative true abundance and abundance indices for each year and simulation of the projection standardized to the average abundance or abundance index over time

### PACKAGES ####
library(tidyverse)
library(here)
library(sdmTMB)
source(here("R", "sim_stratmean_fn.R"))



### DATA SET UP ####
# Directories
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
plots <- here("outputs", "plots")

# Parameters
species <- "summerflounder"
season  <- "fall"
ages      <- 0:7
years     <- 1:15
nsims   <- 1:100
ids     <- sprintf("%03d", nsims)

#Data
pop <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_abund-dist.rds", species, season, .x))))
survdat_sq <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_sq_survey.rds", species, season, .x))))
survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))
survdat_reall <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_reall_survey.rds", species, season, .x))))
dist          <- map(ids, ~readRDS(here(dist.dat, sprintf("%s_%s_%s_dist-only.rds", species, season, .x))))




# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

# find the total survey area
survey_area <- as.integer(sum(strata_wts$Area_SqNm))


#function to set selectivity
source(here("R/selectivity_fns.R"))
q = force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1)
(selectivity_values <- q(ages))


## True Abundance ####
# calculate the relative true abundance from the simulated population and distribution
trueN <- map(pop, ~as_tibble(.$N) |>
                mutate(age = ages) |>
                pivot_longer(cols = all_of(years),
                             names_to = "year",
                             values_to = "N") |>
                mutate(N = N * selectivity_values[as.character(age)]) |>     # apply selectivity to obtain surveyed pop
                summarise(N = sum(N), .by = "year") |> # calculate the sum of N across ages
                mutate(rel_N = N/mean(N), # standardize the annual population by the average population size over the projection
                       year = as.integer(year),
                       scenario = "True")
              ) |>
  map_dfr(~pluck(.), .id = "pop")



trueN_yr <- trueN |>
  group_by(year) |>   # Group data by year to calculate yearly summaries across populations
  summarise(
    mean_N = mean(N), # Compute mean abundance (N) across all pops for each year
    sd_N = sd(N),     # Compute sd of abund across pop for each year
    n_pops = n(),     # Count pops
    se_N = sd_N / sqrt(n_pops), # std error of the mean abundance
    cv_N = sd_N / mean_N,   # cv
    log_se_N = sqrt(log(1 + cv_N^2)), # approximate the std error (log space)
    log_lower = log(mean_N) - 1.96 * log_se_N, # Lower 95% confidence limit in log-space
    log_upper = log(mean_N) + 1.96 * log_se_N, # Upper 95% confidence limit in log-space
    ci_lower = exp(log_lower), # Transform lower limit back to original scale (exponentiate)
    ci_upper = exp(log_upper)) # Transform upper limit back to original scale (exponentiate)


ggplot(trueN_yr, aes(x = year, y = mean_N)) +
  geom_line(color = "black", linewidth = 1) +  # Line for mean abundance
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "lightgray", alpha = 0.5) +  # Confidence interval
  labs(title = "True abundance", x = "Year", y = "True Abundance (N)") +
  theme_minimal(base_size = 14)




# calculate the median relative abundance value and the upper and lower confidence intervals across simulations
# true_med <- trueN |>
#   group_by(year, scenario) |>
#   summarise(med = median(rel_N),
#             lower = quantile(rel_N, 0.025),
#             upper = quantile(rel_N, 0.975)) |>
#   mutate(type = "Relative True Abundance")


## Abundance Index ####
# calculate the abundance index and relative abundance index for each of the scenarios

# extract the strata that were used to predict spatial distributions in sdmTMB
strat <- map(dist, ~unique(.$strat))

### Status Quo ####
# ihat_sq <- map(survdat_sq, ~as_tibble(.$setdet) |> filter(strat %in% strat) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
#               mutate(rel_ihat = stratmu/mean(stratmu),
#                      scenario = "Status Quo")) |>
#   map_dfr(~pluck(.), .id = "sim")

#FOR ONE REALIZATION OF THE POPULATION
# ihat_sq1 <- survdat_sq[[1]]$setdet |> as_tibble() |>
#   filter(strat %in% strat) |>
#   group_by(sim, year, strat) |>
#   summarise(towct = length(unique(set)),
#             mu = sum(n)/towct,
#             var = ifelse(towct == 1, 0,
#                          sum((n - mu)^2)/(towct - 1))) |>
#   left_join(strata_wts, by = "strat") |>
#   mutate(wt_mu = Area_SqNm * mu,
#          wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
#   ungroup() |>
#   group_by(sim, year) |>
#   summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
#             stratvar = sum(wt_var),
#             cv = sqrt(stratvar)/stratmu) |>
#   mutate(rel_ihat = stratmu/mean(stratmu),
#          scenario = "Status Quo")

#FOR ALL THE REALIZATIONS OF THE POPULATION
ihat_sq_all <- map2_dfr(survdat_sq, seq_along(survdat_sq), function(surv, pop_num) {
  surv$setdet |>
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(towct = length(unique(set)),
              mu = sum(n)/towct,
              var = ifelse(towct == 1, 0,
                           sum((n - mu)^2)/(towct - 1)),
              .groups = "drop") |>
    left_join(strata_wts, by = "strat") |>
    mutate(wt_mu = Area_SqNm * mu,
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
    group_by(sim, year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area,
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu,
              .groups = "drop") |>
    mutate(scenario = "Status Quo",
           pop = pop_num)})

# group by pop and sims to compute rel_ihat within each population. If it get computed in the previous step, all populations get combined so if I compute rel_ihat as stratmu / mean(stratmu) across all rows, it end up dividing by a different (larger) mean value.
ihat_sq_all <- ihat_sq_all |>
  group_by(pop, sim) |>
  mutate(rel_ihat = stratmu / mean(stratmu),
         log_se = sqrt(log(1 + cv^2)),
         log_lower = log(stratmu) - 1.96 * log_se,
         log_upper = log(stratmu) + 1.96 * log_se,
         ci_lower = exp(log_lower),
         ci_upper = exp(log_upper)) |>
  ungroup()


### Precluded Survey ####
### #FOR ONE REALIZATION OF THE POPULATION
# ihat_precl <- map(survdat_precl, ~as_tibble(.) |> filter(strat %in% strat) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
#                  mutate(rel_ihat = stratmu/mean(stratmu),
#                         scenario = "Preclusion")) |>
#   map_dfr(~pluck(.), .id = "sim")

# ihat_precl1 <- survdat_precl[[1]] |> as_tibble() |>
#   filter(strat %in% strat) |>
#   group_by(sim, year, strat) |>
#   summarise(towct = length(unique(set)),
#             mu = sum(n)/towct,
#             var = ifelse(towct == 1, 0,
#                          sum((n - mu)^2)/(towct - 1))) |>
#   left_join(strata_wts, by = "strat") |>
#   mutate(wt_mu = Area_SqNm * mu,
#          wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
#   ungroup() |>
#   group_by(sim, year) |>
#   summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
#             stratvar = sum(wt_var),
#             cv = sqrt(stratvar)/stratmu) |>
#   mutate(rel_ihat = stratmu/mean(stratmu),
#          scenario = "Preclusion")

#FOR ALL THE REALIZATIONS OF THE POPULATION
ihat_precl_all <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  surv |> #it is already a data.table
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(towct = length(unique(set)),
              mu = sum(n)/towct,
              var = ifelse(towct == 1, 0,
                           sum((n - mu)^2)/(towct - 1)),
              .groups = "drop") |>
    left_join(strata_wts, by = "strat") |>
    mutate(wt_mu = Area_SqNm * mu,
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
    group_by(sim, year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area,
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu,
              .groups = "drop") |>
    mutate(scenario = "Preclusion",
           pop = pop_num)
})

# Now group by pop (and sim if needed) to compute rel_ihat within each population
ihat_precl_all <- ihat_precl_all %>%
  group_by(pop, sim) %>%
  mutate(rel_ihat = stratmu / mean(stratmu),
         log_se = sqrt(log(1 + cv^2)),
         log_lower = log(stratmu) - 1.96 * log_se,
         log_upper = log(stratmu) + 1.96 * log_se,
         ci_lower = exp(log_lower),
         ci_upper = exp(log_upper)) %>%
  ungroup()



### Reallocated Survey ####
# ihat_reall <- map(survdat_reall, ~as_tibble(.x) |> filter(strat %in% strat) |> sim_stratmean(strata_wts = strata_wts, survey_area = survey_area) |>
#                  mutate(rel_ihat = stratmu/mean(stratmu),
#                         scenario = "Reallocation")) |>
#   map_dfr(~pluck(.), .id = "sim")
#

#FOR ONE REALIZATION OF THE POPULATION
# ihat_reall1 <- survdat_reall[[1]] |> as_tibble() |>
#   filter(strat %in% strat) |>
#   group_by(sim, year, strat) |>
#   summarise(towct = length(unique(set)),
#             mu = sum(n)/towct,
#             var = ifelse(towct == 1, 0,
#                          sum((n - mu)^2)/(towct - 1))) |>
#   left_join(strata_wts, by = "strat") |>
#   mutate(wt_mu = Area_SqNm * mu,
#          wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
#   ungroup() |>
#   group_by(sim, year) |>
#   summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
#             stratvar = sum(wt_var),
#             cv = sqrt(stratvar)/stratmu) |>
#   mutate(rel_ihat = stratmu/mean(stratmu),
#          scenario = "Reallocation")

#FOR ALL THE REALIZATIONS OF THE POPULATION
ihat_reall_all <- map2_dfr(survdat_reall, seq_along(survdat_reall), function(surv, pop_num) {
  surv |> #it is already a data.table
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    summarise(towct = length(unique(set)),
              mu = sum(n)/towct,
              var = ifelse(towct == 1, 0,
                           sum((n - mu)^2)/(towct - 1)),
              .groups = "drop") |>
    left_join(strata_wts, by = "strat") |>
    mutate(wt_mu = Area_SqNm * mu,
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) |>
    group_by(sim, year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area,
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu,
              .groups = "drop") |>
    mutate(scenario = "Reallocation",
           pop = pop_num)})

# Now group by pop (and sim if needed) to compute rel_ihat within each population
ihat_reall_all <- ihat_reall_all |>
  group_by(pop, sim) |>
  mutate(rel_ihat = stratmu / mean(stratmu),
         log_se = sqrt(log(1 + cv^2)),
         log_lower = log(stratmu) - 1.96 * log_se,
         log_upper = log(stratmu) + 1.96 * log_se,
         ci_lower = exp(log_lower),
         ci_upper = exp(log_upper)) |>
  ungroup()


### Bind Indices ####
# bind all three scenario dataframes for efficient plotting and statistic calculation
#indices25 <- bind_rows(ihat_sq1, ihat_precl1, ihat_reall1)
indices <- bind_rows(ihat_sq_all,ihat_precl_all,ihat_reall_all)


#indices_old<-readRDS(here(surv.prods, str_c(species, season, "all-ihat-25survs.rds", sep = "_")))
#scup_fall_all-ihat_1pop-25survs


# calculate the median relative abundance indices and the upper and lower confidence intervals across simulations and scenarios
# indices_med <- indices |>
#   group_by(year, scenario) |>
#   summarise(med = median(rel_ihat),
#             lower = quantile(rel_ihat, 0.025),
#             upper = quantile(rel_ihat, 0.975)) |>
#   mutate(type = "Estimated Relative Abundance Index")


## Plots ####

pops_to_plot <- c(2,5,88)
sim_to_plot <- c(1,100)

trueN_sub     <- filter(trueN, pop %in% pops_to_plot)
indices_sub <- filter(indices, pop %in% pops_to_plot)
indices_sub_sim <- filter(indices, pop %in% pops_to_plot, sim %in% sim_to_plot)
trueN_yr_pop_sub <- filter(trueN_yr_pop, pop %in% pops_to_plot)



RelAbundPlot <- ggplot() +
  aes(x = year) +
  geom_line(data = trueN_sub, aes(y = rel_N, color = scenario), linewidth = 1) +
  geom_point(data = indices_sub, aes(y = rel_ihat, color = fct_inorder(scenario)),
             position = position_dodge(width = 0.5)) +
  labs(x = "Year",
       y = "Relative Abundance", title = "Relative abundance of scup", subtitle = "Fall") +
  facet_wrap(~ pop, scales = "free_y", ncol = 1) +
  ylim(0, NA) +
  theme(legend.position = "bottom", legend.title =element_blank(),
        text = element_text(size = 13),  legend.text = element_text(size = 11))



ggsave(str_c(species, season, "RelAbundPlot.png", sep = "_"),
       plot = RelAbundPlot,
       device = "png",
       # last_plot(),
       here(plots),
       width = 8, height = 6)

## SAVE THE DATA ####
saveRDS(trueN, here(surv.prods, str_c(species, season, "rel-TrueN-100pops.rds", sep = "_")))
saveRDS(indices, here(surv.prods, str_c(species, season, "all-ihat-25survs-100pops.rds", sep = "_")))

saveRDS(ihat_sq_all, here(surv.prods, str_c(species, season, "100pops-25sims-sq_rel-ihat.rds", sep = "_")))
saveRDS(ihat_precl_all, here(surv.prods, str_c(species, season, "100pops-25sims-precl_rel-ihat.rds", sep = "_")))
saveRDS(ihat_reall_all, here(surv.prods, str_c(species, season, "100pops-25sims-reall_rel-ihat.rds", sep = "_")))



ggplot() +
  aes(x = year) +
  geom_line(data = trueN_sub, aes(y = rel_N, color = scenario), linewidth = 1) +
  geom_point(data = indices_sub,
             aes(y = rel_ihat, color = fct_inorder(scenario)),
             position = position_dodge(width = 0.5),
             alpha = 0.4, size = 1.5) +
  labs(x = "Year", y = "Relative Abundance",
    title = "Relative Abundance of Scup", subtitle = "Fall") +
  facet_wrap(~ pop, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Status Quo" = "salmon",
                                "Preclusion" = "goldenrod",
                                "Reallocation" = "steelblue",
                                "True" = "grey40")) +
  scale_x_continuous(breaks = 1:15) +
  ylim(0, NA) +
  theme_minimal() +
  theme(legend.position = "right",
    legend.title = element_blank(),
    text = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 5)))







ggplot(trueN, aes(x = year, y = rel_N, group = pop)) +
  geom_line(alpha = 0.3, color = "steelblue") +
  labs(
    title = "Simulated True Relative Abundance",
    subtitle = "Across 100 populations over 15 years",
    x = "Year", y = "Relative Abundance"
  ) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none")





# Summarize true abundance - median
true_summary_med <- trueN %>%
  group_by(year) %>%
  summarise(median_rel_N = median(rel_N), .groups = "drop")

# Summarize true abundance - median
true_summary_avg <- trueN %>%
  group_by(year) %>%
  summarise(mean_rel_N = mean(rel_N), .groups = "drop")

# Summarize estimates (mean or median) for each scenario - median
index_summary_med <- indices %>%
  group_by(year, scenario) %>%
  summarise(median_rel_ihat = median(rel_ihat), .groups = "drop")

# Summarize estimates (mean or median) for each scenario - mean
index_summary_avg <- indices %>%
  group_by(year, scenario) %>%
  summarise(mean_rel_ihat = mean(rel_ihat), .groups = "drop")

# Plot all
ggplot() +
  geom_line(data = true_summary_med, aes(x = year, y = median_rel_N), color = "black", size = 1.2, linetype = "dashed") +
  geom_line(data = index_summary_med, aes(x = year, y = median_rel_ihat, color = scenario), size = 1) +
  labs(title = "Median Relative Abundance Over Time",
       subtitle = "Dashed line = Simulated true abundance | Colored lines = Estimated abundance",
       y = "Relative Abundance",
       x = "Year") +
  theme(legend.position = "bottom", text = element_text(size = 13))


# Plot all
ggplot() +
  geom_line(data = true_summary_avg, aes(x = year, y = mean_rel_N), color = "black", size = 1.2, linetype = "dashed") +
  geom_line(data = index_summary_avg, aes(x = year, y = mean_rel_ihat, color = scenario), size = 1) +
  labs(title = "Mean Relative Abundance Over Time",
       subtitle = "Dashed line = Simulated true abundance | Colored lines = Estimated abundance",
       y = "Relative Abundance",
       x = "Year") +
  theme(legend.position = "bottom", text = element_text(size = 13))


# Reorder scenarios before plotting
indices <- indices %>%
  mutate(scenario = factor(scenario, levels = c("Status Quo", "Preclusion", "Reallocation")))

# Plot boxplots + dashed line
ggplot() +
  geom_boxplot(data = indices,
               aes(x = factor(year), y = rel_ihat, fill = scenario),
               outlier.size = 0.5, width = 0.7, alpha = 0.8,
               position = position_dodge(width = 0.9)) +
    geom_line(data = true_summary_avg,
            aes(x = year, y = mean_rel_N),
            color = "black", size = 1) +
    scale_fill_manual(values = c(
    "Status Quo"   = "salmon",
    "Preclusion"   = "#DAA520",
    "Reallocation" = "#4682B4"
  )) +
    labs(title = "Estimated Relative Abundance vs. Simulated Truth",
       subtitle = "Solid line = Simulated true abundance | Boxes = Survey-based estimates",
       x = "Year",
       y = "Relative Abundance") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size = 13)) +
  scale_x_discrete(breaks = as.character(1:15))  # show all years







# Summarize by year and scenario
ihat_summary <- indices %>%
  group_by(year, scenario) %>%
  summarise(
    mean_rel_ihat = mean(rel_ihat, na.rm = TRUE),
    sd_rel_ihat = sd(rel_ihat, na.rm = TRUE),
    n = n(),
    se_rel_ihat = sd_rel_ihat / sqrt(n),
    cv_rel_ihat = sd_rel_ihat / mean_rel_ihat,
    log_se = sqrt(log(1 + cv_rel_ihat^2)),
    log_lower = log(mean_rel_ihat) - 1.96 * log_se,
    log_upper = log(mean_rel_ihat) + 1.96 * log_se,
    ci_lower = exp(log_lower),
    ci_upper = exp(log_upper),
    .groups = "drop" # don't forget this if using summarise
  )



ggplot(ihat_summary, aes(x = year, y = mean_rel_ihat, color = scenario, fill = scenario)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, color = NA) +
  labs(title = "Estimated Relative Abundance",
       y = "Relative Abundance",
       x = "Year") +
  theme_minimal() +
  scale_color_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue"
  )) +
  scale_fill_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue"
  )) +
  theme(legend.position = "right",
        text = element_text(size = 14))



# Faceted plot by scenario
ggplot(ihat_summary, aes(x = year, y = mean_rel_ihat, color = scenario, fill = scenario)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.3, color = NA) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_color_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue"
  )) +
  scale_fill_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue"
  )) +
  labs(title = "Estimated Relative Abundance per Survey Scenario",
       subtitle = "",
       y = "Relative Abundance",
       x = "Year") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )




# Calculate mean and uncertainty of RELATIVE TRUE ABUNDANCE
trueN_yr_rel <- trueN %>%
  group_by(year) %>%
  summarise(
    mean_rel_N = mean(rel_N),        # mean relative N
    sd_rel_N   = sd(rel_N),           # standard deviation of relative N
    n_pops     = n(),                 # number of pops (should be 100)
    se_rel_N   = sd_rel_N / sqrt(n_pops),   # SE for relative N
    cv_rel_N   = sd_rel_N / mean_rel_N,     # CV for relative N
    log_se_rel_N = sqrt(log(1 + cv_rel_N^2)), # Delta method SE (log scale)
    log_lower = log(mean_rel_N) - 1.96 * log_se_rel_N,
    log_upper = log(mean_rel_N) + 1.96 * log_se_rel_N,
    ci_lower = exp(log_lower),        # lower 95% CI back to normal scale
    ci_upper = exp(log_upper)         # upper 95% CI back to normal scale
  ) |>
  mutate(scenario = "True")


library(dplyr)

# Make sure columns match (rename if necessary)
trueN_yr_rel <- trueN_yr_rel %>%
  rename(rel_ihat = mean_rel_N)  # Rename for consistency

# Bind rows
combined_data <- bind_rows(ihat_summary, trueN_yr_rel)


ggplot(combined_data, aes(x = year, y = rel_ihat, color = scenario, fill = scenario)) +
  # Confidence ribbons (only where ci_lower and ci_upper exist)
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              data = combined_data,
              inherit.aes = TRUE,
              alpha = 0.3, color = NA) +
    # Lines
  geom_line(size = 1.2) +

  labs(
    title = "True vs Estimated Relative Abundance",
    x = "Year",
    y = "Relative Abundance",
    fill = "Scenario",
    color = "Scenario"
  ) +

  scale_color_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue",
    "True" = "black"
  )) +

  scale_fill_manual(values = c(
    "Status Quo" = "salmon",
    "Preclusion" = "goldenrod",
    "Reallocation" = "steelblue",
    "True" = "black"
  )) +

  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )




