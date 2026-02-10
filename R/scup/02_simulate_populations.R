### created: 01/18/2024
### updated: 04/11/2025

# 02 - SIMULATE POPULATIONS ####


## Objective ####
# Simulate multiple realizations of a single species according to respective life history parameters and recent stock assessment estimates of recruitment and fishing mortality. vlall the populations under the base case scenario for each of the key species of interest, and save the R objects for calling into sim_survey script.

# Outputs: simulated abundance for a user provided species of interested over a set of iterations

### PACKAGES ####
suppressPackageStartupMessages(library(tidyverse))
library(SimSurvey)
library(here)
library(plotly)
source(here("R", "sim_pop_fn.R"))
set.seed(131) # random seed number for generating different populations but for reproducibility of simulation



### DATA SET UP ####

# location of where to save data
pop.dat <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")

#read species data

species <- "scup" # name of species to be simulated
nsims <- 1:100 # number of simulations of the population
ages <- 0:7 # ages to simulate
years <- 1:15 # number of years of a given population


#Recruitment
Rec <- c(107,142,75,61,112)*1e6 # Average Rec from benchmark SA
R_mean <- mean(Rec)
R_sd <- sd(Rec)  # Standard deviation of recruitment
R_cv <- R_sd / R_mean  # Coefficient of Variation

# Derive sigma_R for log-normal variability
sigma_R <- sqrt(log(1 + R_cv^2))  # Conversion from CV to log-scale

set.seed(42)  # For reproducibility
# Simulate recruitment with log-normal variability
Rec_age0 <- R_mean * exp(rnorm(years, mean = 0, sd = sigma_R))

# Fixed fishing mortality
F_fixed <- c(0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.05, 0.05)  # F for ages 0-7
F <- matrix(F_fixed, nrow = length(ages), ncol = length(years)) # Constant F across years
M <- 0.2 # natural mortality
Z <- F + M #total mortality

# Von Bertalanffy growth parameters for both male and female
Linf = 46.6
K = 0.15


# generate a set random numbers to set seed with each iteration for reproducibility
set.seed(42)
seed <- sample.int(1e6, length(nsims))


## SIMULATE ABUNDANCE ####
pop <- sim_pop(iter = length(nsims), ages, years, Rec_age0, Z, Linf, K)

### EXTRACT N AT AGE MATRIX ####
Nage <- map_df(pop, ~pluck(., "N") |>
 # pop[[1]]$N |>
  as_tibble() |>
  mutate(age = ages) |>
  pivot_longer(cols=all_of(years),
               names_to = "year",
               values_to = "Nage"),
               .id = "sim") |>
  janitor::clean_names() |>
  rename(Nage = nage) |>
  mutate(sim = as.integer(sim),
         year = as.integer(year))



#Plot 1
ggplot(Nage, aes(x = year, y = age, fill = Nage)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D", name = "Numbers at Age") + # Color scale
  labs(title = "Numbers at Age",
       x = "Year", y = "Age") +
  #facet_wrap(~sim, ncol = 2) + # Facet by simulation
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom", legend.text = element_text(size = 10, angle = 45))


#Plot 2 - recreates surface plot
nage_data <- Nage  # Rename to avoid conflicts

#Creates a surface plot for each simulation
plots <- nage_data |>
  group_by(sim) |>
  group_split() |>
  lapply(function(nage_data) {
    sim_id <- unique(nage_data$sim) # Get simulation ID
    # Reshape data to a grid format
    surface_data <- nage_data |>
      select(year, age, Nage) |>
      pivot_wider(names_from = year, values_from = Nage, names_prefix = "Year_") |>
      column_to_rownames(var = "age")
    z_matrix <- as.matrix(surface_data)  # Convert to matrix for surface plot
    plot_ly(
      x = unique(nage_data$year), y = unique(nage_data$age), z = z_matrix,
      type = "surface", colorscale = "Viridis") |>
      layout(
        title = paste("Simulation", sim_id),
        scene = list(
          xaxis = list(title = "Year"),
          yaxis = list(title = "Age"),
          zaxis = list(title = "Numbers at Age")))
  }
  )

# Display a specific plot (e.g., first simulation)
plots[[1]]
plots[[2]]
plots[[100]]


library(htmlwidgets)

for (i in seq_along(plots)) {
  saveWidget(
    plots[[i]],
    file = sprintf("simulation_%03d_surface.html", i),
    selfcontained = TRUE
  )
}


# Save each pop object individually
for (i in seq_along(pop)) {
  saveRDS(
    pop[[i]],
    file = file.path(pop.dat, sprintf("%s_pop_%03d.rds", species, i))
  )
  message(sprintf("Saved %s_pop_%03d.rds", species, i))
}


# Save each Nage object individually
sim_ids <- unique(Nage$sim)

# Loop through and save each one
for (i in sim_ids) {
  nage_i <- Nage |> filter(sim == i)
  saveRDS(nage_i, file = file.path(Nage.dat, sprintf("scup_Nage_%03d.rds", i)))
  message(sprintf("Saved scup_Nage_%03d.rds", i))
}

## SAVE THE DATA ####
saveRDS(pop, here(pop.dat, str_c(species, length(nsims), "pop100.rds", sep = "_")))
saveRDS(Nage, here(Nage.dat, str_c(species, length(nsims), "Nage100.rds", sep = "_")))
