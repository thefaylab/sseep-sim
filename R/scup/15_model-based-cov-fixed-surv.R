# ============================================================
# SCENARIO 6
# Random strat in WEAs
# Supplemental fixed stations as a covariate
# ============================================================

#Packages
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

# DATA SET UP
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
mods.data <- here("data", "rds", "surv-prods", "mods_data", "scup")
plots <- here("outputs", "plots")

species <- "scup"
season <- "fall"
ages <- 0:7
years <- 1:15
nsims <- 1:100
ids <- sprintf("%03d", nsims)


# LOAD DATA
precl   <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))

supp <- map(ids, function(id) {
  x <- readRDS(here(survdat, sprintf("%s_%s_%s_25_supp_fixed_survey.rds",
                                     species, season, id)))
  out <- x$setdet #load only setdet data
  rm(x); gc()
  out
})

# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)


precl <- map(precl, ~ .x |> mutate(supplemental_fixed = 0, division = 1))
supp  <- map(supp, ~ .x |> mutate(supplemental_fixed = 1))


com_data <- map2(precl, supp, ~ bind_rows(.x, .y) |>
                   arrange(sim, year, strat, set))

com_data[[1]]
# Add lon/lat
add_latlon <- function(df, crs_proj = 32618) {
  df$x_m <- df$x * 1000
  df$y_m <- df$y * 1000
  sf_obj <- st_as_sf(df, coords = c("x_m", "y_m"), crs = crs_proj, remove = FALSE)
  sf_obj <- st_transform(sf_obj, crs = 4326)
  coords <- st_coordinates(sf_obj)
  df$lon <- coords[, 1]
  df$lat <- coords[, 2]
  return(df)
  }


com_data_ll <- lapply(com_data, function(df) {
  df <- as.data.frame(df)
  df <- add_latlon(df)        # append lon/lat
  return(df)
})



tb_com <- map2_dfr(com_data_ll, seq_along(com_data_ll), function(surv, pop_num) {
  surv |>
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    left_join(strata_wts, by = "strat")|>
    mutate(scenario = "Sc6",
           pop = pop_num,
           SEASON = "FALL") |>
    rename(YEAR = year,
           X = x,
           Y = y,
           DECDEG_LAT = lat,
           DECDEG_LON = lon,
           AVGDEPTH = depth,
           STRATUM = strat) |>
    select(set,sim,YEAR,STRATUM,AVGDEPTH,X,Y,DECDEG_LAT,DECDEG_LON,cell,AREA_CODE,N,n,scenario,pop,SEASON, supplemental_fixed)
})



for (p in unique(tb_com$pop)) {
  for (s in unique(tb_com$sim)) {

    sub <- tb_com %>%
      filter(pop == p, sim == s)

    file_name <- sprintf("suppcov_pop%03d_sim%02d.rds", p, s)
    file_path <- file.path(mods.data, file_name)

    write_rds(sub, file_path, compress = "gz")
  }
}


###############################################
##### SCENARIO 6 - SUPPLEMENTAL COVARIATE #####
###############################################
strata <- readRDS(here("data", "rds", "active_strata.rds")) %>%
  select(STRATUM, Region)

mab_strata <- c(3450, 1050, 1610, 1090, 3410, 3380, 3020, 3460 ,3050, 3440, 3260, 3350, 8510,
                1010, 1060, 3080, 3230, 3320, 3290, 8500, 1650, 1690, 7520, 1100, 3110, 3140,
                3170, 1020, 1740, 1700, 1730, 3200, 1660, 1620, 1110, 1070, 1030, 1750, 1710,
                1670, 1630, 8520, 1120, 1080, 1040, 1760, 1720, 1680, 1640, 8530)
length(mab_strata)

sseep.sim.dir <- "C:/Users/croman1/Desktop/UMassD/sseep-sim"
mods.data.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","mods_data", "scup")
fit.out.dir <- file.path(sseep.sim.dir, "data", "rds", "surv-prods","fit_out", "scup")

################ MODEL FIT ####################

pops <- 16:20
sims <- 1:25

jobs <- expand.grid(pop = pops, sim = sims) |>
  arrange(pop, sim)

plan(multisession, workers = 2)

run_fit <- function(i, jobs, mods.data.dir, fit.out.dir, mab_strata) {

  pop <- jobs$pop[i]
  sim <- jobs$sim[i]

  in_path  <- file.path(mods.data.dir, sprintf("suppcov_pop%03d_sim%02d.rds", pop, sim))
  out_path <- file.path(fit.out.dir,  sprintf("sdmTMB_tw_suppcov_pop%03d_sim%02d.rds", pop, sim))

  simdat <- readRDS(in_path)

  simdat_filt <- simdat |>
    filter(AVGDEPTH <= 75, STRATUM %in% mab_strata) |>
    mutate(
      YEAR = as.factor(YEAR),
      AREA_CODE = as.factor(AREA_CODE),
      supplemental_fixed = as.factor(supplemental_fixed)
    )

  TMB::openmp(n = 1)

  simdat_mesh <- make_mesh(
    simdat_filt, c("X", "Y"),
    fmesher_func = fmesher::fm_mesh_2d_inla,
    cutoff = 30,
    max.edge = c(200, 400),
    offset = c(25, 100)
  )

  fit <- sdmTMB(
    n ~ poly(AVGDEPTH, 2) + YEAR + supplemental_fixed - 1,
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

results_suppcov <- future_lapply(
  X = seq_len(nrow(jobs)),
  FUN = run_fit,
  jobs = jobs,
  mods.data.dir = mods.data,
  fit.out.dir = fit.out.dir,
  mab_strata = mab_strata,
  future.seed = TRUE
)

plan(sequential)


