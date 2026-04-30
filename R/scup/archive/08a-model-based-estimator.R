### created: 07/21/2025
### modified: 03/17/2026


# 01 - IMPORT SCUP DATA FROM SURVEY ####

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
mods.data.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","mods_data", "scup")
fit.out.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","fit_out", "scup")


strata <- readRDS(here("data", "rds", "active_strata.rds")) %>%
  select(STRATUM, Region)

mab_strata <- c(3450, 1050, 1610, 1090, 3410, 3380, 3020, 3460 ,3050, 3440, 3260, 3350, 8510,
                1010, 1060, 3080, 3230, 3320, 3290, 8500, 1650, 1690, 7520, 1100, 3110, 3140,
                3170, 1020, 1740, 1700, 1730, 3200, 1660, 1620, 1110, 1070, 1030, 1750, 1710,
                1670, 1630, 8520, 1120, 1080, 1040, 1760, 1720, 1680, 1640, 8530)
length(mab_strata)



########################################
##### MODEL 1 - NO WIND CONVARIATE #####
########################################

################ MODEL FIT #############

pops <- 1:50
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

plan(multisession, workers = 12)

run_fit <- function(i, jobs, mods.data.dir, fit.out.dir, mab_strata) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  in_path  <- file.path(mods.data.dir, sprintf("sq_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(fit.out.dir,  sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))

    simdat <- readRDS(in_path)

    simdat_filt <- simdat |>
      filter(AVGDEPTH <= 75, STRATUM %in% mab_strata) |>
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
  mab_strata = mab_strata,
  future.seed = TRUE
)



plan(sequential)

################ MODEL PRED #############

fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup"

#plug in grid
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
fall_grid <- readRDS(file = here(sdmtmb.dir, "scup", "data", "scup_fall_grid_122024.rds"))
grid_scup <- fall_grid |> dplyr::select(X, Y, mean_2, Cell_Area, AREA_CODE)  |>
  mutate(AVGDEPTH = mean_2)


pops <- 1:100
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

out_dir <- file.path(fit.dir, "mb_index_sq_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = 12)

run_mb_index <- function(i, jobs, fit.dir, out_dir, grid_scup) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

 fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))
 out_path <- file.path(out_dir, sprintf("mb_index_sq_pop%03d_sim%02d.rds", pop, sim))

 fit <- readRDS(fit_path)

 yrs <- unique(fit$data$YEAR)

 pred_grid <- tidyr::crossing(grid_scup, YEAR = factor(yrs, levels = yrs)) |>
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
  grid_scup = grid_scup,
  future.seed = TRUE
)

plan(sequential)



################ SAVE INDEX #############

species <- "scup"
season <- "fall"

mb_dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup/mb_index_sq_check"

files <- list.files(mb_dir,
                    pattern = "^mb_index_sq_pop\\d{3}_sim\\d{2}\\.rds$",
                    full.names = TRUE)

mb_all <- map_dfr(files, readRDS) |>
  dplyr::mutate(YEAR = as.numeric(YEAR)) |>
  dplyr::arrange(pop, sim, YEAR)



mb_all |> mutate(YEAR = as.numeric(as.character(YEAR)), id = interaction(pop, sim, drop = TRUE)) |>
  ggplot(aes(x = YEAR, y = est, group = id)) + geom_line(alpha = 0.08) +
  geom_point(alpha = 0.08, size = 0.4) + labs(x = "Year", y = "Model-based index") +
  theme_bw()


model.est <- here("data", "rds", "surv-prods", "fit_out", "scup")
saveRDS(mb_all, here(model.est, str_c(species, season, "model_based_wind_ihat.rds", sep = "_")))









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
  out_path <- file.path(fit.out.dir,  sprintf("sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds", pop, sim))

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
    n ~ poly(AVGDEPTH, 2) + YEAR + AREA_CODE - 1,
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

fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup/fall"

#plug in grid
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
fall_grid <- readRDS(file = here(sdmtmb.dir, "scup", "data", "scup_fall_grid_122024.rds"))
grid_scup <- fall_grid |> dplyr::select(X, Y, mean_2, Cell_Area, AREA_CODE)  |>
  mutate(AVGDEPTH = mean_2)


pops <- 51:100
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

out_dir <- file.path(fit.dir, "mb_index_wind_sq_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = 12)

run_mb_index <- function(i, jobs, fit.dir, out_dir, grid_scup) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(out_dir, sprintf("mb_index_wind_sq_pop%03d_sim%02d.rds", pop, sim))

  fit <- tryCatch(readRDS(fit_path), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  yrs <- unique(fit$data$YEAR)

  pred_grid <- tidyr::crossing(grid_scup, YEAR = factor(yrs, levels = yrs)) |>
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
  grid_scup = grid_scup,
  future.seed = TRUE
)

plan(sequential)



################ SAVE INDEX #############

species = "scup"
season = "fall"
mb_dir_wind <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup/fall/mb_index_wind_sq_check"

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


model.est <- here("data", "rds", "surv-prods", "fit_out", "scup")
saveRDS(mb_all_wind, here(model.est, str_c(species, season, "model_based_wind_ihat.rds", sep = "_")))





########################################
#####  CHECK - PLOTS - MISC CALCS  #####
########################################

##-------------------------------------------------------------------
###This is to check if there are missing model fits
# pops <- 1:50
# sims <- 1:25
#
# jobs <- expand.grid(pop = pops, sim = sims) |>
#   arrange(pop, sim)
#
# jobs$filename <- sprintf(
#   "sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds",
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


# plan(multisession, workers = 12)
#
# future.apply::future_lapply(
#   seq_len(nrow(jobs_to_run)),
#   function(i) run_fit(
#     i,
#     jobs_to_run,
#     mods.data.dir,
#     fit.out.dir,
#     mab_strata
#   )
# )
##-------------------------------------------------------------------


# #Check
#
# set.seed(123)
#
# test_fits <- jobs |>
#   distinct(pop) |>
#   slice_sample(n = 5) |>
#   mutate(sim = 15)
#
#
# fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup"
#
# fits_list <- purrr::map2(
#   test_fits$pop,
#   test_fits$sim,
#   function(pop, sim) {
#
#     fit_path <- file.path(
#       fit.dir,
#       sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim)
#     )
#
#     fit <- readRDS(fit_path)
#
#     cat("\n============================\n")
#     cat("POP:", pop, " SIM:", sim, "\n")
#     cat("============================\n")
#
#     print(fit)
#
#     return(fit)
#   }
# )
#
#
# purrr::map2_dfr(
#   test_fits$pop,
#   test_fits$sim,
#   function(pop, sim) {
#
#     fit_path <- file.path(
#       fit.dir,
#       sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim)
#     )
#
#     fit <- readRDS(fit_path)
#
#     tibble(
#       pop = pop,
#       sim = sim,
#       pdHess = fit$sd_report$pdHess
#     )
#   }
# )
#
# fits_list[[1]]
# fits_list[[2]]
#
# #plug in grid
# sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
# fall_grid <- readRDS(file = here(sdmtmb.dir, "scup", "data", "scup_fall_grid_122024.rds"))
# grid_scup <- fall_grid |> dplyr::select(X, Y, mean_2, Cell_Area, AREA_CODE)  |>
#   mutate(AVGDEPTH = mean_2)
#
#
#
# yrs <- levels(fit$data$YEAR)
#
# #predict
# pred_grid <- tidyr::crossing(grid_scup, YEAR = factor(yrs, levels = yrs)) |>
#   dplyr::mutate(AREA_CODE = factor(AREA_CODE))
# pred <- predict(fit, newdata = pred_grid, type = "response")
# head(pred)
#
#
#
# pred_out <- pred
#
#
# #turn preds into model-based index
# mb_index_sq <- pred_out |>
#   dplyr::group_by(YEAR) |>
#   dplyr::summarise(
#     Ihat_mb = sum(est * Cell_Area, na.rm = TRUE),
#     .groups = "drop"
#   )
# mb_index_sq
#
# mb_index_sq <- pred_out |>
#   dplyr::group_by(pop, sim, YEAR) |>
#   dplyr::summarise(
#     Ihat_mb = sum(est * Cell_Area, na.rm = TRUE),
#     .groups = "drop"
#   )
#
# plot(as.integer(as.character(mb_index_sq$YEAR)), mb_index_sq$Ihat_mb, type = "l",
#      xlab = "Year", ylab = "Model-based index")
#
#
#
# #--------
# check_params <- function(pop, sim, fit.out.dir) {
#
#   fit_path <- file.path(
#     fit.out.dir,
#     sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim)
#   )
#
#   fit <- readRDS(fit_path)
#
#   coefs <- coef(fit)$cond
#
#   data.frame(
#     pop = pop,
#     sim = sim,
#     max_abs_coef = max(abs(coefs), na.rm = TRUE)
#   )
# }
#
# param_diag <- purrr::map2_dfr(
#   jobs$pop,
#   jobs$sim,
#   check_params,
#   fit.out.dir = fit.out.dir
# )
#
#
#
# # ggplot(mb_all,
# #        aes(x = as.factor(YEAR),
# #            y = Ihat_mb,
# #            fill = factor(pop))) +
# #   geom_boxplot(
# #     position = position_dodge(width = 0.75),
# #     width = 0.6,
# #     outlier.shape = NA
# #   ) +
# #   scale_fill_brewer(palette = "Set3") +
# #   labs(
# #     x = "Year",
# #     y = "Model-based index",
# #     fill = "Population"
# #   ) +
# #   theme_bw() +
# #   theme(
# #     legend.position = "top",
# #     axis.text.x = element_text(angle = 45, hjust = 1)
# #   )
#
#
#
# ggplot(mb_all, aes(x = factor(YEAR), y = est)) +
#   geom_boxplot() +
#   labs(x = "Year", y = "Model-based index (ihat_mb)") +
#   theme_bw()
#
#
#
# summ <- mb_all %>%
#   group_by(pop, YEAR) %>%
#   summarise(
#     med = median(Ihat_mb, na.rm = TRUE),
#     lo  = quantile(Ihat_mb, 0.25, na.rm = TRUE),
#     hi  = quantile(Ihat_mb, 0.75, na.rm = TRUE),
#     .groups = "drop"
#   )
#
# ggplot(summ, aes(x = YEAR, group = pop)) +
#   geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.08) +
#   geom_line(aes(y = med), alpha = 0.25) +
#   labs(x = "Year", y = "Model-based index (ihat_mb)") +
#   theme_bw() +
#   theme(legend.position = "none")
#
# ggplot(mb_all, aes(x = factor(YEAR), y = Ihat_mb)) +
#   geom_boxplot(outlier.shape = NA) +
#   scale_y_log10() +
#   labs(x = "Year", y = "Model-based index (ihat_mb)") +
#   theme_bw()
