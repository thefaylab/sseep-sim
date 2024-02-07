#adapted from sseep-analysis to use output columns from simsurvey

### LOAD DATA ####
#test <- readRDS(here("data", "rds", "merged_data_complete.rds")) |>
  # mutate(EXPCATCHWT = ifelse(is.na(EXPCATCHWT), 0, EXPCATCHWT),
  #        tow_code = str_c(STRATUM, CRUISE6, STATION)) |>
  # filter(SVSPP == 141, EST_YEAR %in% c(2018, 2021))

#strata <- readRDS(here("data", "rds", "strata.rds"))

sim_stratmean <- function(setdet, strata_wts, survey_area){

  # surv_dat <- sim$setdet |> as_tibble()

  individual <- setdet |>
    group_by(strat, year) |>
    summarise(towct = length(unique(set)), # calculate unique tows
              mu = sum(n)/towct, # find the average biomass based on unique tows rather than observations to avoid potential duplication
              var = ifelse(towct == 1, 0, # if the tow count equals 1, then variance about the mean should be 0
                           sum((n - mu)^2)/(towct - 1))) |> # if tow count does not equal 1, then find the variance of biomass
    left_join(strata_wts, by = "strat") |> # add each stratum area and relative weight to the dataset based on STRATUM number
    mutate(wt_mu = Area_SqNm * mu, # part one of the stratified mean formula
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) # part one of the stratified variance formula

  # BTSArea <- as.integer(sum(surv_dat$cell_area))

  stratified <- individual |>
    group_by(year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu)  # part two of the stratified variance formula

  return(stratified)

  #return(individual)
}

#GF 2024-02-07
#new version that keeps calculations internal.
#Note that this just performs one calculation, looping over simulations/years/scenarios/species would all be done outside the function
sim_stratmean2 <- function(setdet, strata_wts){

  # surv_dat <- sim$setdet |> as_tibble()

  #figure out the strata to use for this calculation, will be the strata from the survey data frame
  strat_use <- unique(setdet$strat)

  #adjust the weights for the strata used in the calculation
  strata_wts_use <- strata_wts |>
    filter(strat %in% strat_use) |>
    mutate(RelWt = Area_SqNm / sum(Area_SqNm))
  #total survey area for this calculation
  survey_area <- sum(strata_wts_use$Area_SqNm)

  #now calculate the stratifiend random mean and variance

  individual <- setdet |>
    group_by(strat) |>
    summarise(towct = length(unique(set)), # calculate unique tows
              mu = sum(n)/towct, # find the average biomass based on unique tows rather than observations to avoid potential duplication
              var = ifelse(towct == 1, 0, # if the tow count equals 1, then variance about the mean should be 0
                           sum((n - mu)^2)/(towct - 1))) |> # if tow count does not equal 1, then find the variance of biomass
    left_join(strata_wts_use, by = "strat") |> # add each stratum area and relative weight to the dataset based on STRATUM number
    mutate(wt_mu = Area_SqNm * mu, # part one of the stratified mean formula
           wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) # part one of the stratified variance formula

  # BTSArea <- as.integer(sum(surv_dat$cell_area))

  stratified <- individual |>
    #group_by(year) |>
    summarise(stratmu = (sum(wt_mu)) / survey_area, # part two of the stratified mean formula
              stratvar = sum(wt_var),
              cv = sqrt(stratvar)/stratmu)  # part two of the stratified variance formula

  return(stratified)

  #return(individual)
}

# strata.mean <- function(x) {
#
#   straified <- x |>
#     group_by("STRATUM", "EST_YEAR", "SVSPP") |>
#     summarise(stratmu = (sum(wt_mu)) / BTSArea, # part two of the stratified mean formula
#               stratvar = sum(wt_var))  # part two of the stratified variance formula
#
#  return(stratified)
#
# }


#strata.test <- strata.mean(test)

#read in data to compare against test results
#nowind_stratmu <- readRDS(here("data", "rds", "strat-mu_included.rds"))




