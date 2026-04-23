


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

# Plot mean trajectory with uncertainty across populations
ggplot(trueN_stat, aes(x = year, y = mean_rel_N)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line() +
  theme_bw()


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
  crossing(sim = sort(unique(estimates$sim)))

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
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c(
    "Status Quo"       = "#2F397A",
    "Preclusion"       = "#7391BD",
    "Reallocation"     = "#894846",
    "Ratio estimator"  = "#785838",
    "Model based"      = "#93995C",
    "Model based wind" = "#4F6009"
  )) +
  labs(x = "Year", y = "Coverage", color = "Scenario") +
  theme_bw()



### SANITY CHECK: SCALE COMPARISON ####
# Check that estimates and truth are on comparable scales
cov_dat |>
  summarise(min_est = min(estimate, na.rm = TRUE),
            max_est = max(estimate, na.rm = TRUE),
            min_true = min(true_rel, na.rm = TRUE),
            max_true = max(true_rel, na.rm = TRUE),
            .by = scenario)
