#' Stratified mean calculation
#'
#' Calculates stratified relative abundance indices from simulated survey data.
#'
#' @param surv_list List of survey datasets (one per population realization).
#' @param strata_wts Data frame containing strata areas used for weighting.
#' @param survey_area Total survey area used to calculate stratum weights.
#' @param scenario_name Character string identifying the survey scenario.
#' @param value_col Name of the column containing the observation values to summarize (e.g., catch per tow).
#' @param years Optional vector of years to include in the calculation.
#'
#' @return A data frame containing stratified mean abundance (`stratmu`),
#' variance (`stratvar`), coefficient of variation (`cv`), and identifiers
#' for simulation, year, population, and scenario.

calc_stratmean <- function(surv_list,
                      strata_wts,
                      survey_area,
                      scenario_name,
                      value_col = "n",
                      years = NULL) {

  map2_dfr(surv_list, seq_along(surv_list), function(surv, pop_num) {

    dat <- surv |>
      as_tibble() |>
      filter(strat %in% unlist(strat))

    if (!is.null(years)) {
      dat <- dat |> filter(year %in% years)
    }

    dat |>
      group_by(sim, year, strat) |>
      summarise(
        towct = length(unique(set)),
        mu = sum(.data[[value_col]]) / towct,
        var = ifelse(towct < 2, 0, var(.data[[value_col]])),
        .groups = "drop"
      ) |>
      left_join(strata_wts, by = "strat") |>
      mutate(
        W = Area_SqNm / survey_area,
        wt_mu = W * mu,
        wt_var = ifelse(towct < 2, 0, (W^2) * var / towct)
      ) |>
      group_by(sim, year) |>
      summarise(
        stratmu = sum(wt_mu),
        stratvar = sum(wt_var),
        cv = sqrt(stratvar) / stratmu,
        .groups = "drop"
      ) |>
      mutate(
        scenario = scenario_name,
        pop = pop_num
      )
  })
}
