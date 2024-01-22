# Function of species population simulation based on SimSurvey package
#
#
# @param iter: # of simulation
# @param ages: ages of simulated species
# @param years: population projection years
# @param Rec_age0: recruitment values at age 0, see sim_abundance() helper for additional information
# @param Linf:
# @param K:
#
sim_pop <- function(iter, ages, years, Rec_age0, Z, Linf, K, ...) {

  set.seed(round(rnorm(1, mean = 1000, sd = 120), 0))
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
