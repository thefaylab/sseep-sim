### created: 01/18/2024
### updated: 02/07/2025

# 02 - SIMULATE POPULATIONS ####


## Objective ####
# Simulate multiple realizations of a single species according to respective life history parameters and recent stock assessment estimates of recruitment and fishing mortality. vlall the populations under the base case scenario for each of the key species of interest, and save the R objects for calling into sim_survey script.

# Outputs: simulated abundance for a user provided species of interested over a set of iterations

### PACKAGES ####
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(here)
library(plotly)
source(here("R", "sim_pop_fn.R"))
# random seed number for generating different populations but for reproducibility of simulation
set.seed(131)


### DATA SET UP ####
#read species data


# location of where to save data
pop.dat <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")


### name of species to be simulated
species <- "scup"

### number of simulations
nsims <- 1:2

# generate a set random numbers to set seed with each iteration for reproducibility
seed <- sample.int(1e6, length(nsims))


#### STOCK ASSESSMENT VALUES ####
### replace with the values respective to species being simulated

# ages to simulate
ages <- 0:7

# number of years of a given population
years <- 1:15

# Base recruitment value (constant recruitment assumption)
R_base <- 100e6  # 100 million recruits

# Process variability for recruitment
sigmaR <- 0.427  # Recruitment variability, estimated in SR-scup.R

# Generate recruitment values with variability
Rec_age0 <- R_base * exp(rnorm(length(years), mean = 0, sd = sigmaR))

# Fixed fishing mortality
F_fixed <- c(0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.05, 0.05)  # F for ages 0-7
F <- matrix(F_fixed, nrow = length(ages), ncol = length(years))

# natural mortality
M <- 0.2

#total mortality
Z <- F + M

# Von Bertalanffy growth parameters for both male and female
Linf = 46.6
K = 0.15


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


#plot_surface(pop[1], mat = "N")
#


#Plot 1
ggplot(Nage, aes(x = year, y = age, fill = Nage)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D", name = "Numbers at Age") + # Color scale
  labs(title = "Numbers at Age",
       x = "Year", y = "Age") +
  facet_wrap(~sim, ncol = 2) + # Facet by simulation
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


## SAVE THE DATA ####
saveRDS(pop, here(pop.dat, str_c(species, length(nsims), "pop.rds", sep = "_")))
saveRDS(Nage, here(Nage.dat, str_c(species, length(nsims), "Nage.rds", sep = "_")))
