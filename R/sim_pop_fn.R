#' Simulate multiple realizations of a species population
#'
#'
#' @param seed: # of simulation
#' @param ages: ages of simulated species
#' @param years: population projection years
#' @param Rec_age0: recruitment values at age 0, see sim_abundance() helper for additional information
#' @param Linf:
#' @param K:
#'
sim_pop <- function(iter, ages, years, Rec_age0, Z, Linf, K, seed = sample.int(1e6, 1L), ...) {

  set.seed(seed)
  pop <- map(1:iter, ~sim_abundance(ages = ages,
                       years = years,
                       Z = sim_Z(log_mean = log(Z), # changed this from 0.4 to 0.3
                                 log_sd = 0.001,
                                 plot = FALSE),
                       R = sim_R(log_mean = log(Rec_age0),
                                 log_sd = sd(log(Rec_age0)),
                                 plot = FALSE),
                       N0 = sim_N0(N0 = "exp",
                                   plot = FALSE),
                       growth = sim_vonB(Linf = Linf,
                                         K = K,
                                         plot = FALSE)))
}


#' Sample prediction years from sdmTMB outputs
#'
#'
#' Sample locations used in the model fitting process to serve as potential future locations for simulating response values with `sdmTMB_simulate()`
#'
#' @param nyears number of years to sample
#' @param year_data integer vector of year values to sample
#' @param
#' @param mod_preds output from `sdmTMB::predict()`; estimated spatial random effects in link space at a given time and location
#' @param seed Random number seed
#'
#' @returns
#' An input dataframe for sdmTMB_simulate() containing X and Y columns in UTM coordinates, the fixed effects components of the fitted linear predictor, and the estimated spatial random effects


sample_years <- function(nyears, year_data, replace = TRUE, seed = sample.int(1e6, 1L)){
  # create one data frame of future locations
  set.seed(seed)

  year_samps <- map(nyears, ~sample(unique(year_data), size = 1, replace = replace))

  return(year_samps)

}

#'
#'
#'
#' Filter sdmTMB predicted distributions
#'
#'
#'
#'
#' @param year_samps
#' @param preds
#' @param seed
#' @param
#' @param
#' @param
#'
filter_distributions <- function(#nyears, year_data,
  year_samps, preds, seed = sample.int(1e6, 1L)){

  set.seed(seed)

  preds_list <- list(preds)

  future_dists <- map2(preds_list, year_samps, ~filter(.x, EST_YEAR %in% .y)) |>
    map2_dfr(seq(year_samps), ~mutate(., year = rep(.y, length.out = nrow(.x))))

  return(future_dists)
}
#'
#'

