
#pred check
#

out_dir <- file.path(fit.dir, "mb_index_sq_check")

indices_list <- purrr::map2(
  test_fits$pop,
  test_fits$sim,
  function(pop, sim) {

    file_path <- file.path(
      out_dir,
      sprintf("mb_index_sq_pop%03d_sim%02d.rds", pop, sim)
    )

    idx <- readRDS(file_path)

    cat("\n============================\n")
    cat("POP:", pop, " SIM:", sim, "\n")
    cat("============================\n")

    print(idx)

    return(idx)
  }
)


out_dir <- file.path(fit.dir, "mb_index_sq_check")
files <- list.files(out_dir, pattern = "mb_index_sq_pop.*\\.rds$", full.names = TRUE)
mb_index_all <- map_dfr(files, function(f) {

  df <- readRDS(f)

  # extract pop and sim from filename
  name <- basename(f)

  pop <- as.integer(stringr::str_extract(name, "(?<=pop)\\d{3}"))
  sim <- as.integer(stringr::str_extract(name, "(?<=sim)\\d{2}"))

  df |>
    mutate(pop = pop, sim = sim)
})

mb_index_all |>
  mutate(
    YEAR = as.numeric(as.character(year)),
    id = interaction(pop, sim, drop = TRUE)
  ) |>
  ggplot(aes(x = YEAR, y = est, group = id)) +
  geom_line(alpha = 0.08) +
  geom_point(alpha = 0.08, size = 0.4) +
  labs(x = "Year", y = "Model-based index") +
  theme_bw()


mb_index_all2 <- mb_index_all |>
  rename(year = YEAR) |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(est, na.rm = TRUE),  #   cv = se_natural / est,
    cv = sqrt(exp(se^2) - 1),
    var_mean_ihat = sum(se_natural^2, na.rm = TRUE) / (n_years^2),   # this is the model-based analog of var_mean_ihat

    # covariance between annual estimate and the across-year mean
    cov_est_mean = (se_natural^2) / n_years,
    rel_ihat = est / mean_ihat,
    rel_var =
      (se_natural^2 / (mean_ihat^2)) +
      ((est^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * est * cov_est_mean / (mean_ihat^3)),
    rel_var = pmax(rel_var, 0),
    rel_se = sqrt(rel_var),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) |>
  ungroup()


threshold <- quantile(mb_index_all$est, 0.999, na.rm = TRUE)

outliers <- mb_index_all |>
  filter(est > threshold)


threshold <- 5e14
mb_index_trim <- mb_index_all |>
  filter(est < threshold)
nrow(mb_index_all) - nrow(mb_index_trim)



mb_index_trim2 <- mb_index_trim |>
  rename(year = YEAR) |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(est, na.rm = TRUE),
    rel_ihat = est / mean_ihat
  ) |>
  ungroup()


mb_index_trim2  |> mutate(
  year = as.numeric(year),
  id = interaction(pop, sim)
) |>
  ggplot(aes(x = year, y = rel_ihat, group = id)) +
  geom_line(alpha = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw()



mb_index_all |>
  mutate(
    year = as.numeric(year),
    id = interaction(pop, sim)
  ) |>
  ggplot(aes(x = year, y = rel_ihat, group = id)) +
  geom_line(alpha = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw()


mb_index_all |>
  mutate(YEAR = as.numeric(as.character(YEAR))) |>
  ggplot(aes(x = YEAR, y = est, group = interaction(pop, sim))) +
  geom_line(alpha = 0.05)



library(tidyverse)

pops <- 1:5
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

read_mb_index <- function(pop, sim, out_dir) {

  file_path <- file.path(
    out_dir,
    sprintf("mb_index_sq_pop%03d_sim%02d.rds", pop, sim)
  )

  readRDS(file_path)
}

mb_index_all <- purrr::map2_dfr(
  jobs$pop,
  jobs$sim,
  read_mb_index,
  out_dir = out_dir
)

mb_index_all <- mb_index_all |>
  mutate(YEAR = as.numeric(as.character(YEAR))) |>
  arrange(pop, sim, YEAR)


ggplot(mb_index_all, aes(x = factor(YEAR), y = est)) +
  geom_boxplot() +
  facet_wrap(~pop) +
  labs(x = "Year", y = "Estimated index")

ggplot(mb_index_all, aes(x = factor(YEAR), y = se)) +
  geom_boxplot() +
  labs(x = "Year", y = "SE")


mb_index_all |>
  group_by(YEAR) |>
  summarise(
    mean_est = mean(est, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE)
  ) |>
  ggplot(aes(x = YEAR, y = mean_est)) +
  geom_line() +
  geom_point() +
  labs(x = "Year", y = "Mean estimated index")


mb_index_all |>
  mutate(id = paste(pop, sim, sep = "_")) |>
  ggplot(aes(x = YEAR, y = est, group = id)) +
  geom_line(alpha = 0.1) +
  labs(x = "Year", y = "Estimated index")


ggplot(mb_index_all, aes(x = factor(YEAR), y = est)) +
  geom_boxplot()

ggplot(mb_index_all, aes(x = factor(YEAR), y = se)) +
  geom_boxplot()

ggplot(mb_index_all, aes(x = factor(YEAR), y = upr - lwr)) +
  geom_boxplot()


