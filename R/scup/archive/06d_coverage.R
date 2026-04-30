


### PACKAGES ####
library(tidyverse)
library(here)
library(sdmTMB)
library(SimSurvey)
source(here("R", "sim_stratmean_fn.R"))



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

### TRUE POPULATION SUMMARY (FOR PLOTTING ONLY) ####
trueN_stat <- trueN |>
  group_by(year) |>
  summarise(
    mean_rel_N = mean(rel_N),
    sd_rel_N = sd(rel_N),
    n_pops = n(),
    se_rel_N = sd_rel_N / sqrt(n_pops),
    cv_rel_N = sd_rel_N / mean_rel_N,
    sdlog = sqrt(log(1 + cv_rel_N^2)),
    meanlog = log(mean_rel_N) - 0.5 * sdlog^2,
    ci_lower = qlnorm(0.025, meanlog = meanlog, sdlog = sdlog),
    ci_upper = qlnorm(0.975, meanlog = meanlog, sdlog = sdlog)
  )


### Prepare true FOR COVERAGE ####
true <- trueN |>
  transmute(
    pop = as.integer(pop),  # ensure matching type with estimates
    year = as.integer(year),
    true_rel = rel_N  # true relative abundance (same scale as estimates)
  )


### PREPARE ESTIMATES ####
estimates <- indices |>
  transmute(
    pop = as.integer(pop),
    sim = as.integer(sim), # survey replicate
    year = as.integer(year),
    scenario,                   # scenario label
    estimate = rel_ihat,       # estimated relative index
    lower = rel_ci_lower,
    upper = rel_ci_upper
  )


### EXPAND true ACROSS SIM REPLICATES ####
# True population is identical across sims, so repeat it
true_expanded <- true |>
 crossing(sim = sort(unique(estimates$sim))) #unnecesary because i have the join

# Check dimensions
nrow(true)          # should be 1500 = 100 pops * 15 years
nrow(true_expanded) # should be 37500 = 100 * 15 * 25 sims

### JOIN TRUE AND ESTIMATES ####
cov_dat <- estimates |>
  left_join(
    true_expanded,
    by = c("pop", "sim", "year")
  )


sum(is.na(cov_dat$true_rel)) # Check that join worked (no missing truth values) should be 0


### COMPUTE ####
cov_dat <- cov_dat |> mutate(covered = true_rel >= lower & true_rel <= upper)

### Overall
coverage_overall <- cov_dat |>
  summarise(coverage = mean(covered, na.rm = TRUE),
            n = sum(!is.na(covered)),
            .by = scenario)



### Cov by year
coverage_year <- cov_dat |>
  summarise(coverage = mean(covered, na.rm = TRUE),
            n = sum(!is.na(covered)),
            .by = c(scenario, year))

coverage_year <- coverage_year |>
  mutate(scenario = factor(scenario, levels = c(
      "Status Quo","Preclusion", "Reallocation","Ratio estimator","Model based","Model based wind")))


p_cov <- ggplot(coverage_year, aes(x = year, y = coverage, color = scenario)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("Status Quo"       = "#127088",
                                "Preclusion"       = "#C85729",
                                "Reallocation"     = "#92874B",
                                "Ratio estimator"  = "#CD8A39",
                                "Model based"      = "#AC3414",
                                "Model based wind" = "#57643C")) +
  labs(x = "Year", y = "Coverage", color = "Scenario", title = "") +
  theme(text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x  = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold"))



ggsave(str_c(species, season, "coverage", sep = "_"),
       plot = p_cov,
       device = "png",
       # last_plot(),
       here(plots),
       width = 7, height = 5)

### SANITY CHECK: SCALE COMPARISON ####
# Check that estimates and truth are on comparable scales
cov_dat |>
  summarise(min_est = min(estimate, na.rm = TRUE),
            max_est = max(estimate, na.rm = TRUE),
            min_true = min(true_rel, na.rm = TRUE),
            max_true = max(true_rel, na.rm = TRUE),
            .by = c(scenario, year)) |>
  arrange(scenario, year)

cov_dat |>
  filter(scenario %in% c("Model based", "Model based wind")) |>
  group_by(scenario, year) |>
  summarise(
    coverage = mean(covered, na.rm = TRUE),
    mean_est = mean(estimate, na.rm = TRUE),
    mean_true = mean(true_rel, na.rm = TRUE),
    bias = mean(estimate - true_rel, na.rm = TRUE),
    mean_lower = mean(lower, na.rm = TRUE),
    mean_upper = mean(upper, na.rm = TRUE),
    prop_true_above = mean(true_rel > upper, na.rm = TRUE),
    prop_true_below = mean(true_rel < lower, na.rm = TRUE),
    .groups = "drop"
  )



ihat_model_all2 |>
  group_by(year) |>
  summarise(mean_est = mean(est, na.rm = TRUE))


pops<- mb_all |>
  group_by(pop, sim) |>
  summarise(
    mean_ihat = unique(mean_ihat)
  ) |>
  summarise(
    min_mean = min(mean_ihat),
    max_mean = max(mean_ihat),
    sd_mean = sd(mean_ihat)
  )


ihat_model_all2 |>
  group_by(year) |>
  summarise(
    mean_rel = mean(rel_ihat),
    sd_rel = sd(rel_ihat)
  )


cov_dat |>
  summarise(
    # ESTIMATES
    est_min = min(estimate, na.rm = TRUE),
    est_q05 = quantile(estimate, 0.05, na.rm = TRUE),
    est_q25 = quantile(estimate, 0.25, na.rm = TRUE),
    est_q50 = quantile(estimate, 0.50, na.rm = TRUE),
    est_q75 = quantile(estimate, 0.75, na.rm = TRUE),
    est_q95 = quantile(estimate, 0.95, na.rm = TRUE),
    est_max = max(estimate, na.rm = TRUE),

    # TRUTH
    true_min = min(true_rel, na.rm = TRUE),
    true_q05 = quantile(true_rel, 0.05, na.rm = TRUE),
    true_q25 = quantile(true_rel, 0.25, na.rm = TRUE),
    true_q50 = quantile(true_rel, 0.50, na.rm = TRUE),
    true_q75 = quantile(true_rel, 0.75, na.rm = TRUE),
    true_q95 = quantile(true_rel, 0.95, na.rm = TRUE),
    true_max = max(true_rel, na.rm = TRUE),

    .by = c(scenario, year)
  ) |>
  arrange(scenario, year) |>
  print(n = Inf)



cov_dat |>
  filter(scenario == "Model based") |>
  arrange(desc(estimate)) |>
  select(pop, sim, year, estimate) |>
  slice_head(n = 20)




true_abs <- trueN |>
  transmute(
    pop = as.integer(pop),
    year = as.integer(year),
    true_N = N
  )

mb_abs <- mb_all |>
  transmute(
    pop = as.integer(pop),
    sim = as.integer(sim),
    year = as.integer(YEAR),
    scenario = "Model based",
    est_abs = est
  ) |>
  left_join(true_abs, by = c("pop", "year"))

mb_abs |>
  summarise(
    est_min = min(est_abs, na.rm = TRUE),
    est_q05 = quantile(est_abs, 0.05, na.rm = TRUE),
    est_q25 = quantile(est_abs, 0.25, na.rm = TRUE),
    est_q50 = quantile(est_abs, 0.50, na.rm = TRUE),
    est_q75 = quantile(est_abs, 0.75, na.rm = TRUE),
    est_q95 = quantile(est_abs, 0.95, na.rm = TRUE),
    est_max = max(est_abs, na.rm = TRUE),

    true_min = min(true_N, na.rm = TRUE),
    true_q05 = quantile(true_N, 0.05, na.rm = TRUE),
    true_q25 = quantile(true_N, 0.25, na.rm = TRUE),
    true_q50 = quantile(true_N, 0.50, na.rm = TRUE),
    true_q75 = quantile(true_N, 0.75, na.rm = TRUE),
    true_q95 = quantile(true_N, 0.95, na.rm = TRUE),
    true_max = max(true_N, na.rm = TRUE),

    .by = c(scenario, year)
  ) |>
  arrange(scenario, year) |>
  print(n = Inf)



true_abs <- trueN |>
  transmute(
    pop = as.integer(pop),
    year = as.integer(year),
    true_N = N
  )

est_abs_all <- indices |>
  transmute(
    pop = as.integer(pop),
    sim = as.integer(sim),
    year = as.integer(year),
    scenario,
    est_abs = ihat   # ← THIS is the absolute index for all scenarios
  ) |>
  left_join(true_abs, by = c("pop", "year"))


est_abs_all |>
  summarise(
    est_min = min(est_abs, na.rm = TRUE),
    est_q05 = quantile(est_abs, 0.05, na.rm = TRUE),
    est_q25 = quantile(est_abs, 0.25, na.rm = TRUE),
    est_q50 = quantile(est_abs, 0.50, na.rm = TRUE),
    est_q75 = quantile(est_abs, 0.75, na.rm = TRUE),
    est_q95 = quantile(est_abs, 0.95, na.rm = TRUE),
    est_max = max(est_abs, na.rm = TRUE),

    true_min = min(true_N, na.rm = TRUE),
    true_q05 = quantile(true_N, 0.05, na.rm = TRUE),
    true_q25 = quantile(true_N, 0.25, na.rm = TRUE),
    true_q50 = quantile(true_N, 0.50, na.rm = TRUE),
    true_q75 = quantile(true_N, 0.75, na.rm = TRUE),
    true_q95 = quantile(true_N, 0.95, na.rm = TRUE),
    true_max = max(true_N, na.rm = TRUE),

    .by = c(scenario, year)
  ) |>
  arrange(scenario, year) |>
  print(n = Inf)
