### created: 07/21/2025
### modified: 02/04/2026


# 01 - IMPORT SCUP DATA FROM SURVEY ####

## OBJECTIVE ####
# prepare data for sdmTMB model fits and predicting


## LOAD PACKAGES ####
suppressPackageStartupMessages(library(tidyverse))
library(here)
library(sf)
library(sdmTMB)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringr)
library(readr)
library(dplyr)
# library(marmap)
# library(raster)

here()

## LOAD DATA ####
sseep.sim.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim"
mods.data.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","mods_data", "scup")
fit.out.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","fit_out", "scup")

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
  geom_point(aes(x = AVGDEPTH, y = N)) +
 # facet_wrap(~YEAR) +
  labs(x = "Depth (m)", y = "Catch (weight) per tow", subtitle = "Fall")


mab_strata <- c(3450, 1050, 1610, 1090, 3410, 3380, 3020, 3460 ,3050, 3440, 3260, 3350, 8510, 1010,
                1060, 3080, 3230, 3320, 3290, 8500, 1650, 1690, 7520, 1100, 3110, 3140, 3170, 1020,1740,1700,1730,3200,
                1660,1620,1110,1070,1030,1750,1710,1670,1630,8520,1120,1080,1040,1760,1720,1680,1640,8530)
length(mab_strata)




### filter specific strata ####

simdat_filt <- simdat_1_1 |> filter(AVGDEPTH <= 75,  STRATUM %in% mab_strata) |>
  mutate(YEAR = as.factor(YEAR)) |>
  group_by(set, sim, SEASON, YEAR)


ggplot(simdat_filt) +
  geom_point(aes(x = AVGDEPTH, y = N)) +
  # facet_wrap(~YEAR) +
  labs(x = "Depth (m)", y = "Catch (weight) per tow", subtitle = "Fall")


simdat_mesh <- make_mesh(simdat_filt, xy_cols = c("X", "Y"), cutoff = 10)
plot(simdat_mesh)



simmod_tw <- sdmTMB(N ~ poly(AVGDEPTH,2) + YEAR - 1,
                    data = simdat_filt,
                    mesh = simdat_mesh,
                    family = tweedie(link = "log"),
                    spatial = "on",
                    time = "YEAR",
                    spatiotemporal = "IID",
                    control = sdmTMBcontrol(newton_loops = 1),
                    silent = FALSE)

sanity(simmod_tw)
tidy(simmod_tw)
tidy(simmod_tw, effects = "ran_pars")







# fall_mod_tw <- readRDS(here("sdmtmb", "scup", "data","mods","comps", "m7d_fall_tw.rds"))
# tidy(fall_mod_tw)
# tidy(fall_mod_tw, effects = "ran_pars")
#
#
#


pops <- 1:10
sims <- 1:25

# Log file (CSV) to track success/failure + timing
log_file <- file.path(fit.out.dir, "fit_log_sq_scup.csv")
if (!file.exists(log_file)) {
  writeLines("pop,sim,status,elapsed_sec,n_rows,message", con = log_file)
}

## Loop for ultiple sims and pops ####
for (pop in pops) {
  for (sim in sims) {

    in_name  <- sprintf("sq_pop%03d_sim%02d.rds", pop, sim)
    in_path  <- file.path(mods.data.dir, in_name)

    out_name <- sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim)
    out_path <- file.path(fit.out.dir, out_name)

    # Skip if already fit
    if (file.exists(out_path)) next

    # If input missing, log and continue
    if (!file.exists(in_path)) {
      cat(sprintf("%d,%d,missing_input,NA,NA,%s\n", pop, sim, shQuote(in_path)),
          file = log_file, append = TRUE)
      next
    }

    start_time <- Sys.time()

    # Wrap everything so one failure doesn't kill the whole run
    res <- tryCatch({

      simdat <- readRDS(in_path)

      simdat_filt <- simdat %>%
        filter(AVGDEPTH <= 75, STRATUM %in% mab_strata) %>%
        mutate(YEAR = as.factor(YEAR)) %>%
        group_by(set, sim, SEASON, YEAR)

      # Basic guard: if filtering leaves nothing, skip
      if (nrow(simdat_filt) == 0) stop("No rows after filtering (AVGDEPTH/STRATUM).")

      simdat_mesh <- make_mesh(simdat_filt, xy_cols = c("X", "Y"), cutoff = 10)

      fit <- sdmTMB(
        N ~ poly(AVGDEPTH, 2) + YEAR - 1,
        data = simdat_filt,
        mesh = simdat_mesh,
        family = tweedie(link = "log"),
        spatial = "on",
        time = "YEAR",
        spatiotemporal = "IID",
        control = sdmTMBcontrol(newton_loops = 1),
        silent = TRUE   # set FALSE if you want console output
      )

      # Save fit object to disk (and ONLY to disk)
      saveRDS(fit, out_path)

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat(sprintf("%d,%d,ok,%.3f,%d,%s\n",
                  pop, sim, elapsed, nrow(simdat_filt), "saved"),
          file = log_file, append = TRUE)

      # IMPORTANT: drop fit from memory explicitly
      rm(fit, simdat, simdat_filt, simdat_mesh)
      gc()

      TRUE

    }, error = function(e) {

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      # Escape commas in message for CSV safety
      msg <- gsub(",", ";", conditionMessage(e))
      cat(sprintf("%d,%d,error,%.3f,NA,%s\n", pop, sim, elapsed, msg),
          file = log_file, append = TRUE)

      rm(list = intersect(c("simdat","simdat_filt","simdat_mesh"), ls()))
      gc()

      FALSE
    })

    # optional progress print
    message(sprintf("pop %03d sim %02d -> %s", pop, sim, if (isTRUE(res)) "OK" else "FAIL"))
  }
}



#Check for 1 pop
fit.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup"

pop <- 1
sim <- 1
fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))

fit <- readRDS(fit_path)

fit$formula

#plug in grid
sdmtmb.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis/sdmtmb"
fall_grid <- readRDS(file = here(sdmtmb.dir, "scup", "data", "scup_fall_grid_122024.rds"))
grid_scup <- fall_grid |> dplyr::select(X, Y, mean_2, Cell_Area)  |>
  mutate(AVGDEPTH = mean_2)

yrs <- levels(fit$data$YEAR)

#predict
pred_grid <- tidyr::crossing(grid_scup, YEAR = factor(yrs, levels = yrs))
pred <- predict(fit, newdata = pred_grid, type = "response")
head(pred)

pred_out <- pred


#turn preds into model-based index
mb_index <- pred_out |>
  dplyr::group_by(YEAR) |>
  dplyr::summarise(
    Ihat_mb = sum(est * Cell_Area, na.rm = TRUE),
    .groups = "drop"
  )
mb_index

plot(as.integer(as.character(mb_index$YEAR)), mb_index$Ihat_mb, type = "l",
     xlab = "Year", ylab = "Model-based index")


#All pops x sims
pops <- 1:10
sims <- 1:25

out_dir <- file.path(fit.dir, "mb_index_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (pop in pops) {
  for (sim in sims) {

    fit_path <- file.path(fit.dir, sprintf("sdmTMB_tw_sq_pop%03d_sim%02d.rds", pop, sim))
    if (!file.exists(fit_path)) next

    fit <- readRDS(fit_path)

    yrs <- levels(fit$data$YEAR)

    pred_grid <- tidyr::crossing(grid_scup, YEAR = factor(yrs, levels = yrs))
    pred_out  <- predict(fit, newdata = pred_grid, type = "response")

    mb_index <- pred_out |>
      group_by(YEAR) |>
      summarise(Ihat_mb = sum(est * Cell_Area, na.rm = TRUE), .groups = "drop") |>
      mutate(pop = pop, sim = sim)

    saveRDS(mb_index, file.path(out_dir, sprintf("mb_index_pop%03d_sim%02d.rds", pop, sim)))

    rm(fit, pred_grid, pred_out, mb_index)
    gc()
  }
}


#Plot

mb_dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim/data/rds/surv-prods/fit_out/scup/mb_index_check"

files <- list.files(mb_dir,
                    pattern = "^mb_index_pop\\d{3}_sim\\d{2}\\.rds$",
                    full.names = TRUE)

mb_all <- map_dfr(files, readRDS)


ggplot(mb_all,
       aes(x = as.factor(YEAR),
           y = Ihat_mb,
           fill = factor(pop))) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6,
    outlier.shape = NA
  ) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Year",
    y = "Model-based index",
    fill = "Population"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

