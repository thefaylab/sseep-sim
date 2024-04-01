#' Simulate a logisitic selectivity curve with forced values at age
#'
#' @description Simulate a logisitic selectivity based on the steepness and sigmoid's midpoint, but allows a given age to be forced to a given value of selectivity
#'
#' @param k: the steepness of the curve; inherited from `sim_logistic()`
#' @param x0: the x-value of the sigmoind's midpoint; inherited from `sim_logistic()`
#' @param plot: plot relationship; inherited from `sim_logistic()`
#' @param force_age: whether a selectivity at age will be forced to a given value; default is FALSE
#' @param age: which age will have a forced selectivity
#' @param force_sel: the value that the selecivity will be given
#'
#' @returns
#' A function to simulate a logistic selectivity for a set of ages that returns a vector containing the selectivities at age
#'
#'
#'@examples
#'
#' q_fun <- force_sim_logistic(k = 1.2, x0 = 2, plot = TRUE, force_age = TRUE, age = 0, force_sel = 0)
#'
#' # with no dimensions
#' sel1 <- q_fun(x = 0:7)
#'
#' # with multiple columns
#' sel2 <- q_fun(replicate(2, 0:7))
#'
#' # with no dimensions but mutiple iterations of an age
#' sel3 <- q_fun(rep(0:7, 2))
#'
#'

force_sim_logistic <- function(k = 2, x0 = 3, plot = FALSE, force_age = FALSE, age, force_sel = 0) {
  function(x = NULL) {
    y <- 1 / (1 + exp(-k * (x - x0)))
    if (plot) plot(x, y, type = "b")
    if (force_age) { # if force_age = TRUE
      if (is.null(dim(y)) == TRUE){ # test if y lacks dimensions
      names(y) <- x # make named number
      #if (force_age){ #
        y[which(names(y) == as.character(age))] <- force_sel # find the instances where the named number is equal to the age provided by user to force the selectivity and replace its value with the force_sel value provided
        y # return
      } else { # if y has more than one column of selectivities
      final_cols <- dim(y)[2] # extract the number of columns in y that will be filled with selectivity values
      colnames(y) <- seq(ncol(y)) # assign column names with the values of the total number columns
      y <- cbind(y, x) # add a column for ages
      #if (force_age){ # if force_age = TRUE
        y[x == age] <- force_sel # find the instances where the x value provided to the function is equal to the age provided by the user to force the selectivity and replace its value with the force_sel value provided
      y[,seq(final_cols)] # only return the selectivities
    }
      } else {
      y
      }
  }
}

