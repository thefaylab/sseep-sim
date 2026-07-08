### created: 01/18/2024
### updated: 07/07/2026


# PIPELINE OVERVIEW:
# BLOCK 1: build a fixed set of tow locations (wind strata + their non-wind
#          siblings) using ONE population (pop 1) as the template design.
# BLOCK 2: apply those SAME fixed locations to every population in a chunk,
#          simulating the supplemental survey without any new random sampling.
# BLOCK 3: combine the standard/preclusion survey with the supplemental
#          fixed-wind survey into one dataset per population.
# BLOCK 4: compute stratmean, then calculate the abundance index.

### PACKAGES ####
library(sdmTMB)
library(SimSurvey)
library(purrr)
library(data.table)
library(here)
library(dplyr)
library(tidyverse)
library(cubelyr)
suppressPackageStartupMessages(library(tidyverse))
source(here("R/selectivity_fns.R"))
set.seed(131)


### DATA SET UP ####
# Directories
pop.dat   <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")
dist.dat  <- here("data", "rds", "dists")
survdat   <- here("data", "rds", "survdat")
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"


# Parameters
species   <- "scup"
season    <- "fall"
ages      <- 0:7
years     <- 1:15
nsims     <- 1:100
nsurveys  <- 25
years_post <- 6:15 # Define the years affected by the survey change
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))


# Trawl survey configuration
trawl_dim       <- c(2.7, 0.014)
catch_q         <- force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1)
set_den         <- 0.001
min_sets        <- 3
age_sampling    <- "stratified"
age_space_group <- "set"
resample_cells  <- TRUE


# BLOCK 1
# only population 1 is used here. This defines a single, fixed
# survey design (tow locations) that gets applied to every other population
# in BLOCK 2. This is not a loop over populations by design.

i <- 1

pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

pop_inside <- pop_i

# 1. Identify wind cells
wind <- pop_i$grid[["AREA_CODE"]] == 1



# 2. Restratify: give wind cells their own unique stratum code
# sim_survey() has no direct argument for controlling where it samples. Tt allocates tows (via min_sets/set_den) based entirely on the "strat"
# attribute already stored in the grid object.
# Wind cells currently share a stratum code with the surrounding non-wind area, so sampling inside the wind footprint is left to chance.
# Splitting wind cells into their own stratum here makes min_sets/set_den apply directly to them, guaranteeing tow coverage inside the wind area.


strat_layer <- pop_i$grid[["strat"]]
offset <- 10000 ## offset must exceed the max existing stratum code (here, max = 3450) to guarantee new wind-substrata codes never collide with real stratum codes

new_strat <- strat_layer
wind_idx <- which(as.vector(wind))
new_strat[wind_idx] <- as.vector(strat_layer)[as.vector(wind_idx)] + offset

pop_inside$grid[["strat"]] <- new_strat


# expect these two counts to be EQUAL - confirms NAs in `wind` and
# `strat_layer` line up on the same out-of-domain cells (land/no-data),
# not a mismatch introduced by this step
sum(is.na(as.vector(wind)))
sum(is.na(as.vector(strat_layer)))

# sanity check: confirms restratification worked as intended
# diagonal entries (row == column) = non-wind cells that keep their orifinal stratum code (unchanged)
# off-diagonal entries where column = row + offset = wind cells that were reassigned to a new stratum code, split out of their original stratum
# counts should match so if stratum 1010 had 192 wind cells, they now appear all under new stratum (1010 + offset), with the  remaining non-wind cells still counted under 1010

table(as.vector(strat_layer), as.vector(new_strat), useNA = "ifany")

survey_inside_i <- sim_survey(
  pop_inside,
  n_sims = nsurveys,
  trawl_dim = trawl_dim,
  q = catch_q,
  set_den = set_den,
  min_sets = min_sets,
  age_sampling = age_sampling,
  age_space_group = age_space_group,
  resample_cells = TRUE)



# 3. Identify original strata that contain any wind cells needed so the supplemental survey pool includes the FULL stratum (wind + non-wind cells sharing that original stratum), not just the wind cells alone

strata_with_wind <- unique(as.vector(strat_layer)[as.vector(wind)])
strata_with_wind <- strata_with_wind[!is.na(strata_with_wind)]

fixed_locs_wind_strata <- survey_inside_i$setdet |>
  filter(year == 1, sim %in% 1:nsurveys) |>
  filter(strat %in% strata_with_wind | strat %in% (strata_with_wind + offset)) |>  #tows are kept in both the wind strata and the wind + offset strata
  distinct(sim, set, strat, x, y, cell, depth, AREA_CODE,
           tow_area, cell_area, strat_cells, strat_area,
           strat_sets, cell_sets) |>
  mutate(division = 1)

fixed_inside_wind_strata <- fixed_locs_wind_strata |>
  tidyr::crossing(year = years_post) |>
  arrange(sim, year, strat, set) |>
  mutate(set = row_number()) |>   # important
  select(sim, year, strat, x, y, cell, depth, AREA_CODE, division,
         tow_area, cell_area, strat_cells, strat_area,
         strat_sets, cell_sets, set)

saveRDS(fixed_inside_wind_strata, here(survdat, "scup_fall_pop001_sims_year01_fixed_inside_wind_strata.rds"))




fixed_inside_wind_strata |>
  filter(AREA_CODE == 2) |>   # non-wind remainder
  distinct(strat, strat_sets, strat_area)

fixed_inside_wind_strata |>
  filter(AREA_CODE == 1) |>   # wind cells
  distinct(strat, strat_sets, strat_area)

#BLOCK 2
#RUN survey with fixed locss

# Define chunking
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))
chunk_id <- 5 #change from 1 to 5
this_chunk <- chunks[[chunk_id]]


# c) Load fixed tow locations (created here in this script)
#fixed_inside_wind_strata <- readRDS(here("data", "rds", "survdat", "inside_fixed_locs_wind_strata_survey.rds"))
fixed_inside_wind_strata <- readRDS(here("data", "rds", "survdat", "scup_fall_pop001_sims_year01_fixed_inside_wind_strata.rds"))


# d) Loop over populations in this chunk
set.seed(123)
for (i in this_chunk) {
  message(sprintf("Running fixed supp survey inside wind areas for population %03d of chunk %d...", i, chunk_id))

  # Load population abundance-distribution object
  pop_i <- readRDS(here(dist.dat, sprintf("%s_%s_%03d_abund-dist.rds", species, season, i)))

  # Run survey for fixed stations only (no random sampling, no reallocation)
  survey_fixed_wind_strata <- sim_survey(
    sim              = pop_i,
    n_sims           = nsurveys,           # 25 identical supplemental survey replicates
    trawl_dim        = trawl_dim,
    q                = catch_q,
    resample_cells   = FALSE,              # fixed locations only
    custom_sets      = fixed_inside_wind_strata,       # use inside-wind grid
    age_sampling     = "random",
    age_space_group  = "set")

  # Save result
  saveRDS(survey_fixed_wind_strata, here(survdat, sprintf("%s_%s_%03d_%d_supp_fixed_survey_wind_strata.rds", species, season, i, nsurveys)))
  message(sprintf("Saved fixed supp survey for population %03d.", i))
}



# ---- VALIDATION: confirms BLOCK 2 sampled correctly ----
# 1) no tows leaked into strata unrelated to wind areas
survey_fixed_wind_strata$setdet |>
  distinct(strat) |>
  filter(!(strat %in% strata_with_wind | strat %in% (strata_with_wind + offset))) #returns 0. any stratum showing up here would mean a tow leaked in from a stratum that does not have wind cells


# 2) both wind (AREA_CODE 1) and non-wind (AREA_CODE 2) cells present,
#    as intended - the design deliberately samples full strata, not just
#    wind cells

survey_fixed_wind_strata$setdet %>% count(AREA_CODE) #shows area code 1 and 2 but thats ok because the code intentionally was designed for that.


# 3) tow allocation per wind stratum matches expected set_den/min_sets math
survey_fixed_wind_strata$setdet |>
  filter(AREA_CODE == 1) |>
  count(strat) #23 strata (new strata) that are wind affected with their respective amount of tows per sim/year/pop
#750 / 25 / 10 = 3 exactly the min_sets - the majority of the strata
#For 11650 (1500 rows): 1500 / 25 / 10 = 6 tows per sim.
#For 11690 (1250 rows): 1250 / 25 / 10 = 5 tows per sim.


# fixed_inside_wind_strata %>%
#   filter(AREA_CODE == 1) %>%
#   distinct(strat, strat_sets, strat_area)

# 4) only years_post (6-15) appear, confirming years are restricted correctly
table(survey_fixed_wind_strata$setdet$year) #confirms years 6 to 15 only



### BLOCK 3: Arrange supp fixed survey + preclusion ####
#survdat_supp <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_supp_fixed_survey_wind_strata.rds", species, season, .x))))

chunk_size <- 100
chunks <- split(nsims, ceiling(nsims / chunk_size))
chunk_id <- 1
this_chunk <- chunks[[chunk_id]]

for (i in this_chunk) {
  message(sprintf("Combining preclusion + fixed supplemental survey inside wind strata for population %03d...", i))

  precl_survey <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_precl_survey.rds", species, season, i, nsurveys)))
  supp_wa_survey  <- readRDS(here(survdat, sprintf("%s_%s_%03d_%d_supp_fixed_survey_wind_strata.rds", species, season, i, nsurveys)))

  precl_setdet <- precl_survey |>
    mutate(survey_type = "standard_precl")

  # supp_wa_survey is filtered to years_post because the supplemental design
  # only applies after the survey change; precl_survey already only contains
  # the years it's meant to (no filter needed there)
  supp_wa_setdet <- supp_wa_survey$setdet |>
    filter(year %in% years_post) |>
    mutate(survey_type = "supplemental_fixed_wind_area")

  survey_combined <- bind_rows(precl_setdet, supp_wa_setdet) |>
    arrange(sim, year, strat, set)

  saveRDS(survey_combined, here(survdat, sprintf("%s_%s_%03d_%d_supp_wa+precl_survey.rds", species, season, i, nsurveys)))

  message(sprintf("Saved supp_wa + precl survey for population %03d: %d rows", i, nrow(survey_combined)))
}




#### BLOCK 4: Calculate abundance index ####

strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)

survey_area <- as.integer(sum(strata_wts$Area_SqNm))

## Abundance Index ####
source(here("R/stratmean_fn.R"))
sseep.sim <- "D:/UMassD/sseep-sim"
dist.dat  <- here(sseep.sim, "data", "rds", "dists")
survdat   <- here(sseep.sim, "data", "rds", "survdat")

species <- "scup"
season  <- "fall"
ages    <- 0:7
years   <- 1:15
nsims   <- 1:100
ids     <- sprintf("%03d", nsims)

survdat_supp <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_supp_wa+precl_survey.rds", species, season, .x))))


### Supplemental Survey ####

# The offset codes (e.g. 11010) exist only to force sim_survey() to guarantee tow coverage inside wind cells during sampling (BLOCK 1/2). They don't represent real separate strata for estimation - every offset-coded tow still physically belongs to its original survey stratum. So before computing the stratified mean, recode strat back to the original code.

survdat_supp <- map(survdat_supp, ~
                      as_tibble(.x) |>
                      mutate(strat = ifelse(strat > offset, strat - offset, strat)))

# years 1:5 predate the wind supplement
survdat_supp_y1_5 <- map(survdat_supp, ~
                           as_tibble(.x) |>
                           filter(year %in% 1:5))

stratmean_supp_y1_5 <- calc_stratmean(
  surv_list     = survdat_supp_y1_5,
  strata_wts    = strata_wts,
  survey_area   = survey_area,
  scenario_name = "Supplemental Fixed",
  value_col     = "n",
  years         = 1:5)

# years 6:15
survdat_supp_y6_15 <- map(survdat_supp, ~
                            as_tibble(.x) |>
                            filter(year %in% 6:15))

stratmean_supp_y6_15 <- calc_stratmean(
  surv_list     = survdat_supp_y6_15,
  strata_wts    = strata_wts,
  survey_area   = survey_area,
  scenario_name = "Supplemental Fixed",
  value_col     = "n",
  years         = 6:15)

# ---- merge ----
stratmean_supp_all <- bind_rows(stratmean_supp_y1_5, stratmean_supp_y6_15) |>
  arrange(pop, sim, year)

## Abundance Index ####
# calculate the abundance index and relative abundance index for each scenario
# Group to compute rel_ihat
ihat_supp_all <- stratmean_supp_all |>
  group_by(pop, sim) |>
  mutate(
    n_years = n_distinct(year),
    mean_ihat = mean(stratmu, na.rm = TRUE),
    var_mean_ihat = sum(stratvar, na.rm = TRUE) / (n_years^2),
    cov_stratmu_mean = stratvar / n_years,
    rel_ihat = stratmu / mean_ihat,
    rel_var =
      (stratvar / (mean_ihat^2)) +
      ((stratmu^2) * var_mean_ihat / (mean_ihat^4)) -
      (2 * stratmu * cov_stratmu_mean / (mean_ihat^3)),
    rel_se = sqrt(pmax(rel_var,0)),
    rel_cv = rel_se / rel_ihat,
    rel_log_sd = sqrt(log(1 + rel_cv^2)),
    rel_log_mean = log(rel_ihat) - 0.5 * rel_log_sd^2,
    rel_ci_lower = qlnorm(0.025, meanlog = rel_log_mean, sdlog = rel_log_sd),
    rel_ci_upper = qlnorm(0.975, meanlog = rel_log_mean, sdlog = rel_log_sd)
  ) %>%
  ungroup()

saveRDS(ihat_supp_all, here(surv.prods, str_c(species, season, "100pops-25sims-supp_rel-ihat.rds", sep = "_")))

