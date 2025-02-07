#Stock Recruitment
#BH
#
#Table 43. Summary assessment results; Spawning Stock Biomass (SSB) in metric tons (mt);
#Recruitment (R) at age 0 in millions; Fishing Mortality (F) for fully recruited ages 2-7+
#Terceiro, Mark (2012). Stock Assessment of scup (Stenotomus chrysops) for 2012.

ssb <- c(11479, 15031, 14341, 11320, 8602, 7459, 10361, 8413, 6949, 5563,
         4202, 3624, 5412, 5438, 6592, 13340, 27792, 53561, 80358, 104409,
         110325, 120631, 130122, 142113, 163555, 178334, 208869, 209171,
         205496, 199034, 182915)

rec <- c(132, 127, 82, 63, 118, 67, 100, 89, 36, 37,
       61, 35, 29, 78, 97, 222, 146, 138, 84, 84,
       127, 197, 222, 218, 185, 98, 107, 142, 75, 61, 112)

# Define the Beverton-Holt function
beverton_holt <- function(ssb, alpha, beta) {
  return((alpha * ssb) / (1 + beta * ssb))
}

# Fit the model using Nonlinear Least Squares (nls)
bh_fit <- nls(rec ~ beverton_holt(ssb, alpha, beta),
              start = list(alpha = 5, beta = 0.00004),
              trace = TRUE)

# Get estimated parameters
(alpha_est <- coef(bh_fit)["alpha"])
(beta_est <- coef(bh_fit)["beta"])


#Calculate predicted rec
rec_pred <- (alpha_est * ssb) / (1 + beta_est * ssb)

#log scale residuals
residuals <- (log(rec) - log(rec_pred))

#Estimate sigma R
(sigmaR_est <- sd(residuals))  # Standard deviation of recruitment residuals


# Compare actual vs. predicted recruitment
plot(ssb, rec, col = "blue", pch = 16, xlab = "Spawning Stock Biomass (SSB)", ylab = "Recruitment")
points(ssb, rec_pred, col = "red", pch = 16)
legend("topright", legend = c("Observed", "Predicted"), col = c("blue", "red"), pch = 16)

