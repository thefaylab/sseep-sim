### created: 07/21/2025
### modified: 04/24/2026


# 01 - IMPORT SUMMERFLOUNDER DATA FROM SURVEY ####

## OBJECTIVE ####
# sdmTMB model fits and predicting


## LOAD PACKAGES ####
library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringr)
library(readr)
library(dplyr)
library(profvis)
library(doParallel)
library(foreach)
library(future)
library(future.apply)
library(fmesher)
library(TMB)
# library(marmap)
# library(raster)

here()

## LOAD DATA ####
sseep.sim.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim"
mods.data.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","mods_data", "sumflounder", "spring")
fit.out.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","fit_out", "sumflounder","spring")


strata <- readRDS(here("data", "rds", "active_strata.rds")) %>%
  select(STRATUM, Region)


#spring summerflounder
sf_strata <- c(1050, 1130, 3020, 1100, 1060, 1140, 1010, 3140, 1020,
               1070, 1110, 1120, 1730, 1030, 1040, 1740, 1750, 1760,
               1690, 1700, 1710, 1720, 1650, 1660, 1670, 1680, 3380,
               1610, 1620, 1630)


length(sf_strata)

########################################
##### MODEL 1 - NO WIND CONVARIATE #####
########################################

################ MODEL FIT #############

pops <- 1:50 #run in batch of 50
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

plan(multisession, workers = 12)

run_fit <- function(i, jobs, mods.data.dir, fit.out.dir, sf_strata) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  in_path  <- file.path(mods.data.dir, sprintf("sq_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(fit.out.dir,  sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))

  simdat <- readRDS(in_path)

  simdat_filt <- simdat |>
    filter(AVGDEPTH <= 75, STRATUM %in% sf_strata) |>
    mutate(YEAR = as.factor(YEAR), AREA_CODE = as.factor(AREA_CODE))

  TMB::openmp(n = 1)

  simdat_mesh <- make_mesh(
    simdat_filt, c("X", "Y"),
    fmesher_func = fmesher::fm_mesh_2d_inla,
    cutoff = 30,
    max.edge = c(200, 400),
    offset = c(25, 100))

  fit <- sdmTMB(
    n ~ poly(AVGDEPTH, 2) + YEAR  - 1,
    data = simdat_filt,
    mesh = simdat_mesh,
    family = tweedie(link = "log"),
    spatial = "on",
    time = "YEAR",
    spatiotemporal = "IID",
    control = sdmTMBcontrol(newton_loops = 1),
    silent = TRUE
  )

  saveRDS(fit, out_path)

  rm(fit, simdat, simdat_filt, simdat_mesh)
  gc()
}

results_sq <- future_lapply(
  X = seq_len(nrow(jobs)),
  FUN = run_fit,
  jobs = jobs,
  mods.data.dir = mods.data.dir,
  fit.out.dir = fit.out.dir,
  sf_strata = sf_strata,
  future.seed = TRUE
)



plan(sequential)

################ MODEL PRED #############

fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/sumflounder/spring"

#plug in grid
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
spring_grid <- readRDS(file = here(sdmtmb.dir, "sumflounder", "data", "sumf_spring_grid_062025.rds"))
grid_sf <- spring_grid |> dplyr::select(X, Y, mean_2, Cell_Area, AREA_CODE)  |>
  mutate(AVGDEPTH = mean_2)


pops <- 1:100
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

out_dir <- file.path(fit.dir, "mb_index_sq_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = 12)

run_mb_index <- function(i, jobs, fit.dir, out_dir, grid_sf) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(out_dir, sprintf("mb_index_sq_pop%03d_sim%02d.rds", pop, sim))

  fit <- tryCatch(readRDS(fit_path), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  yrs <- unique(fit$data$YEAR)

  pred_grid <- tidyr::crossing(grid_sf, YEAR = factor(yrs, levels = yrs)) |>
    dplyr::mutate(AREA_CODE = factor(AREA_CODE, levels = levels(fit$data$AREA_CODE)))

  pred_out <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)

  mb_index_sq <- get_index(pred_out, area = pred_grid$Cell_Area, bias_correct = TRUE) |>
    dplyr::mutate(pop = pop, sim = sim)

  saveRDS(mb_index_sq, out_path)

  rm(fit, pred_grid, pred_out, mb_index_sq)
  gc()
}


results_mb <- future_lapply(
  X = seq_len(nrow(jobs)),
  FUN = run_mb_index,
  jobs = jobs,
  fit.dir = fit.dir,
  out_dir = out_dir,
  grid_sf = grid_sf,
  future.seed = TRUE
)

plan(sequential)


################ SAVE INDEX #############

species = "summerflounder"
season = "spring"
mb_dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/sumflounder/spring/mb_index_sq_check"

files <- list.files(mb_dir,
                    pattern = "^mb_index_sq_pop\\d{3}_sim\\d{2}\\.rds$",
                    full.names = TRUE)


mb_all <- map_dfr(files, readRDS) |>
  dplyr::mutate(YEAR = as.numeric(as.character(YEAR))) |>
  dplyr::arrange(pop, sim, YEAR)


#PLOT THE INDEX
mb_all |> mutate(YEAR = as.numeric(as.character(YEAR)), id = interaction(pop, sim, drop = TRUE)) |>
  ggplot(aes(x = YEAR, y = est, group = id)) + geom_line(alpha = 0.08) +
  geom_point(alpha = 0.08, size = 0.4) + labs(x = "Year", y = "Model-based index") +
  theme_bw()


mb_all <- mb_all |>
  dplyr::filter(!(pop == 24 & sim == 12)) #in this particular case, pop 24 and sim 12 was the series that looked off so i removed it

#PLOT THE INDEX (again)
mb_all |> mutate(YEAR = as.numeric(as.character(YEAR)), id = interaction(pop, sim, drop = TRUE)) |>
  ggplot(aes(x = YEAR, y = est, group = id)) + geom_line(alpha = 0.08) +
  geom_point(alpha = 0.08, size = 0.4) + labs(x = "Year", y = "Model-based index") +
  theme_bw()

model.est <- here("data", "rds", "surv-prods", "fit_out", "sumflounder")
saveRDS(mb_all, here(model.est, str_c(species, season, "model_based_ihat.rds", sep = "_")))





########################################
#####  MODEL 2  -  WIND CONVARIATE #####
########################################

pops <- 1:50 #run in batch of 50
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

plan(multisession, workers = 12)

run_fit <- function(i, jobs, mods.data.dir, fit.out.dir, sf_strata) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  in_path  <- file.path(mods.data.dir, sprintf("sq_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(fit.out.dir,  sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))

  simdat <- readRDS(in_path)

  simdat_filt <- simdat |>
    filter(AVGDEPTH <= 75, STRATUM %in% sf_strata) |>
    mutate(YEAR = as.factor(YEAR), AREA_CODE = as.factor(AREA_CODE))

  TMB::openmp(n = 1)

  simdat_mesh <- make_mesh(
    simdat_filt, c("X", "Y"),
    fmesher_func = fmesher::fm_mesh_2d_inla,
    cutoff = 30,
    max.edge = c(200, 400),
    offset = c(25, 100))

  fit <- sdmTMB(
    n ~ poly(AVGDEPTH, 2) + YEAR  - 1,
    data = simdat_filt,
    mesh = simdat_mesh,
    family = tweedie(link = "log"),
    spatial = "on",
    time = "YEAR",
    spatiotemporal = "IID",
    control = sdmTMBcontrol(newton_loops = 1),
    silent = TRUE
  )

  saveRDS(fit, out_path)

  rm(fit, simdat, simdat_filt, simdat_mesh)
  gc()
}

results_sq <- future_lapply(
  X = seq_len(nrow(jobs)),
  FUN = run_fit,
  jobs = jobs,
  mods.data.dir = mods.data.dir,
  fit.out.dir = fit.out.dir,
  sf_strata = sf_strata,
  future.seed = TRUE
)

plan(sequential)


################ MODEL PRED #############

fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/sumflounder/spring"

#plug in grid
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
spring_grid <- readRDS(file = here(sdmtmb.dir, "sumflounder", "data", "sumf_spring_grid_062025.rds"))
grid_sf <- spring_grid |> dplyr::select(X, Y, mean_2, Cell_Area, AREA_CODE)  |>
  mutate(AVGDEPTH = mean_2)


pops <- 51:100
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

out_dir <- file.path(fit.dir, "mb_index_wind_sq_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = 12)

run_mb_index <- function(i, jobs, fit.dir, out_dir, grid_sf) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(out_dir, sprintf("mb_index_wind_sq_pop%03d_sim%02d.rds", pop, sim))

  fit <- tryCatch(readRDS(fit_path), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  yrs <- unique(fit$data$YEAR)

  pred_grid <- tidyr::crossing(grid_sf, YEAR = factor(yrs, levels = yrs)) |>
    dplyr::mutate(AREA_CODE = factor(AREA_CODE, levels = levels(fit$data$AREA_CODE)))

  pred_out <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)

  mb_index_sq <- get_index(pred_out, area = pred_grid$Cell_Area, bias_correct = TRUE) |>
    dplyr::mutate(pop = pop, sim = sim)

  saveRDS(mb_index_sq, out_path)

  rm(fit, pred_grid, pred_out, mb_index_sq)
  gc()
}


results_mb <- future_lapply(
  X = seq_len(nrow(jobs)),
  FUN = run_mb_index,
  jobs = jobs,
  fit.dir = fit.dir,
  out_dir = out_dir,
  grid_sf = grid_sf,
  future.seed = TRUE
)

plan(sequential)



################ SAVE INDEX #############

species = "summerflounder"
season = "spring"
mb_dir_wind <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/sumflounder/spring/mb_index_wind_sq_check"

files_wind <- list.files(mb_dir_wind,
                    pattern = "^mb_index_wind_sq_pop\\d{3}_sim\\d{2}\\.rds$",
                    full.names = TRUE)


mb_all_wind <- map_dfr(files_wind, readRDS) |>
  dplyr::mutate(YEAR = as.numeric(as.character(YEAR))) |>
  dplyr::arrange(pop, sim, YEAR)


#PLOT THE INDEX
mb_all_wind |> mutate(YEAR = as.numeric(as.character(YEAR)), id = interaction(pop, sim, drop = TRUE)) |>
  ggplot(aes(x = YEAR, y = est, group = id)) + geom_line(alpha = 0.08) +
  geom_point(alpha = 0.08, size = 0.4) + labs(x = "Year", y = "Model-based index") +
  theme_bw()


model.est <- here("data", "rds", "surv-prods", "fit_out", "sumflounder")
saveRDS(mb_all_wind, here(model.est, str_c(species, season, "model_based_wind_ihat.rds", sep = "_")))





########################################
#####  CHECK - PLOTS - MISC CALCS  #####
########################################

##-------------------------------------------------------------------
# mb_all <- mb_all %>%
#   group_by(pop, sim) %>%
#   mutate(
#     mean_ihat = mean(est),
#     rel_ihat_mb = est / mean_ihat
#   ) %>%
#   ungroup()
#
#
#
###This is to check if there are missing model fits
# pops <- 1:50
# sims <- 1:25
#
# jobs <- expand.grid(pop = pops, sim = sims) |>
#   arrange(pop, sim)
#
# jobs$filename <- sprintf(
#   "sdmTMB_tw_sq_pop%03d_sim%02d.rds",
#   jobs$pop,
#   jobs$sim
# )
#
# jobs$filepath <- file.path(fit.out.dir, jobs$filename)
#
# # only detect missing files
# jobs$needs_run <- !file.exists(jobs$filepath)
#
# jobs_to_run <- jobs[jobs$needs_run, ]
#
# nrow(jobs_to_run)
# head(jobs_to_run)
#
# jobs_to_run$filename
#
#
# plan(multisession, workers = 12)
#
# future.apply::future_lapply(
#   seq_len(nrow(jobs_to_run)),
#   function(i) run_fit(
#     i,
#     jobs_to_run,
#     mods.data.dir,
#     fit.out.dir,
#     sf_strata
#   )
# )

##-------------------------------------------------------------------


# ggplot(mb_all,
#        aes(x = as.factor(YEAR),
#            y = Ihat_mb,
#            fill = factor(pop))) +
#   geom_boxplot(
#     position = position_dodge(width = 0.75),
#     width = 0.6,
#     outlier.shape = NA
#   ) +
#   scale_fill_brewer(palette = "Set3") +
#   labs(
#     x = "Year",
#     y = "Model-based index",
#     fill = "Population"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = "top",
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )



ggplot(mb_all, aes(x = factor(YEAR), y = est)) +
  geom_boxplot() +
  labs(x = "Year", y = "Model-based index (ihat_mb)") +
  theme_bw()



summ <- mb_all %>%
  group_by(pop, YEAR) %>%
  summarise(
    med = median(Ihat_mb, na.rm = TRUE),
    lo  = quantile(Ihat_mb, 0.25, na.rm = TRUE),
    hi  = quantile(Ihat_mb, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(summ, aes(x = YEAR, group = pop)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08) +
  geom_line(aes(y = med), alpha = 0.25) +
  labs(x = "Year", y = "Model-based index (ihat_mb)") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(mb_all, aes(x = factor(YEAR), y = Ihat_mb)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_log10() +
  labs(x = "Year", y = "Model-based index (ihat_mb)") +
  theme_bw()



#--------------
pops <- 1:50
sims <- 1:25

jobs_all <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

out_dir <- file.path(fit.dir, "mb_index_sq_check")

existing_files <- list.files(
  out_dir,
  pattern = "\\.rds$",
  full.names = FALSE
)

done <- tibble(file = existing_files) |>
  mutate(
    pop = as.integer(str_extract(file, "(?<=pop)\\d{3}")),
    sim = as.integer(str_extract(file, "(?<=sim)\\d{2}"))
  ) |>
  filter(!is.na(pop), !is.na(sim)) |>
  distinct(pop, sim)

jobs_missing <- jobs_all |>
  anti_join(done, by = c("pop", "sim"))

nrow(jobs_missing)
jobs_missing

plan(multisession, workers = 12)

results_missing <- future_lapply(
  X = seq_len(nrow(jobs_missing)),
  FUN = run_mb_index,
  jobs = jobs_missing,
  fit.dir = fit.dir,
  out_dir = out_dir,
  grid_sf = grid_sf,
  future.seed = TRUE
)




mb_all |>
  mutate( YEAR = as.numeric(as.character(YEAR)), id = interaction(pop, sim, drop = TRUE)) |>
  ggplot(aes(x = YEAR, y = est, group = id)) +
  geom_line(data = ~ filter(.x, pop != 24), color = "grey70", alpha = 0.08) +
  geom_line(data = ~ filter(.x, pop == 24), color = "red", linewidth = 0.8, alpha = 0.9) +
  labs(x = "Year", y = "Model-based index") + theme_bw()



mb_all |>
  mutate(YEAR = as.numeric(as.character(YEAR)), id = interaction(pop, sim)) |>
  ggplot(aes(YEAR, est, group = id)) +
  geom_line(data = ~ dplyr::filter(.x, !(pop == 24 & sim == 12)), color = "grey50", alpha = 0.1) +
  geom_line(data = ~ dplyr::filter(.x, pop == 24 & sim == 12), color = "red", linewidth = 1.2) +
  theme_bw()

