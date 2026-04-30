################
library(sdmTMB)
library(tidyverse)

set.seed(8675309)

#fit simple cod example from sdmTMB
m <- sdmTMB(density ~ 0 + as.factor(year),
            data = pcod_2011, mesh = pcod_mesh_2011, family = tweedie(link = "log"),
            time = "year",
            #spatial = "off",
            #spatiotemporal = "ar1",
)
# predict
qcs_grid_2011 <- replicate_df(qcs_grid, "year", unique(pcod_2011$year))

p2 <- predict(m, newdata = qcs_grid_2011, return_tmb_object = TRUE,
              bias_correct = TRUE)

# check CV calculations of relative index by generating from the estimated lognormal parameters
# (assumes independence among years)
nsim <- 1000
simulate_check <- get_index(p2) |>
  as_tibble() |>
  mutate(sim = map2(log_est, se, rlnorm, n = nsim)) |>
  unnest(cols = c(sim)) |>
  mutate(rep = rep(1:nsim, n_distinct(year))) |>
  group_by(rep) |>
  mutate(ihat = sim/mean(sim)) |>
  ungroup() |>
  group_by(year) |>
  summarize(median = median(ihat),
            cv = sd(ihat)/mean(ihat))
print(simulate_check)


########################################
# Alternative using get_index_sims()
# generate index using draws from the precision matrix
# this will be slower but technically more accurate because it doesn't assume independence
p <- predict(m, newdata = qcs_grid_2011, nsim = 1000)
x_sims <- get_index_sims(p, return_sims = TRUE)

# obtain CV of relative index using simulation-based approach from get_index_sims()
sim_result <- x_sims |>
  janitor::clean_names() |>
  group_by(iteration) |>
  mutate(ihat = value/mean(value)) |>
  ungroup() |>
  group_by(year) |>
  summarize(median = median(ihat),
            lower = quantile(ihat, 0.025),
            upper = quantile(ihat, 0.975),
            cv = sd(ihat)/mean(ihat))
sim_result

