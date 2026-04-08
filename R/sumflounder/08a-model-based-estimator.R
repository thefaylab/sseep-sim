### created: 07/21/2025
### modified: 04/06/2026


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
mods.data.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","mods_data", "sumflounder")
fit.out.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","fit_out", "sumflounder")

# Example: Load population 1, sim 1
pop <- 1
sim <- 1

file_name <- sprintf("sq_pop%03d_sim%02d.rds", pop, sim)
file_path <- file.path(mods.data.dir, file_name)

# Read the file
simdat_1_1 <- readRDS(file_path)

strata <- readRDS(here("data", "rds", "active_strata.rds")) %>%
  select(STRATUM, Region)


#Catch
ggplot(simdat_1_1) +
  geom_point(aes(x = AVGDEPTH, y = n)) +
 # facet_wrap(~YEAR) +
  labs(x = "Depth (m)", y = "Catch (weight) per tow", subtitle = "Fall")


# mab_strata <- c(3610, 3600, 1200, 3560, 1160, 1190, 1230, 1250, 3450, 3460, 1050, 1090, 1130,
#                    3020, 1100, 3050, 1060, 3080, 1010, 3110, 3140, 3170, 1730, 3200, 3230, 1690,
#                    3260, 3290, 3320, 1650, 3350, 3380, 3410, 1610)

mab_strata <- c(1200, 1160, 1190, 3450, 3460, 1050, 1090, 1130, 3020,
                1100, 3050, 1060, 3080, 1010, 3110, 3140, 3170, 1730,
                3200, 3230, 1690, 3260, 3290, 3320)
length(mab_strata)


pops <- 1:50 #run in batch of 50
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

# Log file (CSV) to track success/failure + timing
log_file <- file.path(fit.out.dir, "fit_log_sq_sumflounder.csv")
if (!file.exists(log_file)) {
  writeLines("pop,sim,status,elapsed_sec,n_rows,message", con = log_file)
}

plan(multisession, workers = 12)

run_fit <- function(i, jobs, mods.data.dir, fit.out.dir, mab_strata) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  in_name  <- sprintf("sq_pop%03d_sim%02d.rds", pop, sim)
  in_path  <- file.path(mods.data.dir, in_name)

  out_name <- sprintf("sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds", pop, sim)
  out_path <- file.path(fit.out.dir, out_name)

  if (!file.exists(in_path)) {
    return(sprintf("%d,%d,missing_input,NA,NA,%s",
                   pop, sim, shQuote(in_path)))
  }

  start_time <- Sys.time()

  res <- tryCatch({

    simdat <- readRDS(in_path)

    simdat_filt <- simdat |>
      filter(AVGDEPTH <= 110, STRATUM %in% mab_strata) |>
      mutate(
        YEAR = as.factor(YEAR),
        AREA_CODE = as.factor(AREA_CODE)
      ) |>
      group_by(set, sim, SEASON, YEAR)

    if (nrow(simdat_filt) == 0) {
      stop("No rows after filtering (AVGDEPTH/STRATUM).")
    }

    TMB::openmp(n = 1)

    simdat_mesh <- make_mesh(
      simdat_filt, c("X", "Y"),
      fmesher_func = fmesher::fm_mesh_2d_inla,
      cutoff = 30,
      max.edge = c(200, 400),
      offset = c(25, 100)
    )

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

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    n_rows <- nrow(simdat_filt)

    out_msg <- sprintf("%d,%d,ok,%.3f,%d,%s",
                       pop, sim, elapsed, n_rows, "saved")

    rm(fit, simdat, simdat_filt, simdat_mesh)
    gc()

    out_msg

  }, error = function(e) {

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    msg <- gsub(",", ";", conditionMessage(e))

    rm(list = intersect(c("simdat","simdat_filt","simdat_mesh"), ls()))
    gc()

    sprintf("%d,%d,error,%.3f,NA,%s",
            pop, sim, elapsed, msg)
  })

  res
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

cat(
  paste(unlist(results_sq), collapse = "\n"),
  file = log_file,
  append = TRUE,
  sep = "\n"
)


plan(sequential)

##-------------------------------------------------------------------
###This is to check if there are missing model fits
pops <- 1:50
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

jobs$filename <- sprintf(
  "sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds",
  jobs$pop,
  jobs$sim
)

jobs$filepath <- file.path(fit.out.dir, jobs$filename)

# only detect missing files
jobs$needs_run <- !file.exists(jobs$filepath)

jobs_to_run <- jobs[jobs$needs_run, ]

nrow(jobs_to_run)
head(jobs_to_run)

jobs_to_run$filename


plan(multisession, workers = 12)

future.apply::future_lapply(
  seq_len(nrow(jobs_to_run)),
  function(i) run_fit(
    i,
    jobs_to_run,
    mods.data.dir,
    fit.out.dir,
    mab_strata
  )
)
##-------------------------------------------------------------------


#Check for 1 pop
fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/sumflounder"

pop <- 1
sim <- 1
fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds", pop, sim))

fit <- readRDS(fit_path)
fit

#plug in grid
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
fall_grid <- readRDS(file = here(sdmtmb.dir, "sumflounder", "data", "sumf_fall_grid_062025.rds"))
grid_sf <- fall_grid |> dplyr::select(X, Y, mean_2, Cell_Area, AREA_CODE)  |>
  mutate(AVGDEPTH = mean_2)

yrs <- levels(fit$data$YEAR)

#predict
pred_grid <- tidyr::crossing(grid_sf, YEAR = factor(yrs, levels = yrs)) |>
  dplyr::mutate(AREA_CODE = factor(AREA_CODE))
pred <- predict(fit, newdata = pred_grid, type = "response")
head(pred)

pred_out <- pred


#turn preds into model-based index
mb_index_sq <- pred_out |>
  dplyr::group_by(YEAR) |>
  dplyr::summarise(
    Ihat_mb = sum(est * Cell_Area, na.rm = TRUE),
    .groups = "drop"
  )
mb_index_sq

mb_index_sq <- pred_out |>
  dplyr::group_by(pop, sim, YEAR) |>
  dplyr::summarise(
    Ihat_mb = sum(est * Cell_Area, na.rm = TRUE),
    .groups = "drop"
  )

plot(as.integer(as.character(mb_index_sq$YEAR)), mb_index_sq$Ihat_mb, type = "l",
     xlab = "Year", ylab = "Model-based index")

#-----------------------------------------------------------------
#Predictions and index - All pops x sims
pops <- 1:100
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

out_dir <- file.path(fit.dir, "mb_index_sq_wind_check")
#dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = 12)

run_mb_index <- function(i, jobs, fit.dir, out_dir, grid_scup) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_wind_sq_pop%03d_sim%02d.rds", pop, sim)) #sq or precl
  if (!file.exists(fit_path)) {
    message(sprintf("pop %03d sim %02d -> MISSING", pop, sim))
    return(NULL)
  }

  out_path <- file.path(out_dir, sprintf("mb_index_sq_wind_pop%03d_sim%02d.rds", pop, sim))

  res <- tryCatch({

    fit <- readRDS(fit_path)

    yrs <- levels(fit$data$YEAR)

    pred_grid <- tidyr::crossing(grid_scup, YEAR = factor(yrs, levels = yrs)) |>
      dplyr::mutate(AREA_CODE = factor(AREA_CODE, levels = levels(fit$data$AREA_CODE)))
    pred_out <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)


    # mb_index_sq <- pred_out |>
    #   dplyr::group_by(YEAR) |>
    #   dplyr::summarise(
    #     Ihat_mb = sum(est * Cell_Area, na.rm = TRUE),
    #     .groups = "drop"
    #   ) |>
    #   dplyr::mutate(pop = pop, sim = sim)
    mb_index_sq <- get_index(pred_out, area = pred_grid$Cell_Area, bias_correct = TRUE) |>
      dplyr::mutate(pop = pop, sim = sim)

    saveRDS(mb_index_sq, out_path)

    rm(fit, pred_grid, pred_out, mb_index_sq)
    gc()

    out_path

  }, error = function(e) {

    rm(list = intersect(c("fit", "pred_grid", "pred_out", "mb_index_sq"), ls()))
    gc()

    NULL
  })

  res
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


#rel index --- change directory depending if its sq or precl

mb_dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup/mb_index_sq_wind_check"

files <- list.files(mb_dir,
                    pattern = "^mb_index_sq_wind_pop\\d{3}_sim\\d{2}\\.rds$",
                    full.names = TRUE)

mb_all <- map_dfr(files, readRDS) |>
  dplyr::mutate(YEAR = as.numeric(YEAR)) |>
  dplyr::arrange(pop, sim, YEAR)



mb_all <- mb_all %>%
  group_by(pop, sim) %>%
  mutate(
    mean_ihat = mean(est),
    rel_ihat_mb = est / mean_ihat
  ) %>%
  ungroup()


model.est <- here("data", "rds", "surv-prods", "fit_out", "scup")
saveRDS(mb_all, here(model.est, str_c(species, season, "model_based_wind_ihat.rds", sep = "_")))


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
