### created: 01/18/2024
### updated: 02/06/2024

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
years <- 1:5

# Recruitment at age 0 for the most recent 5 years, 2015-2019; numbers were reported in millions in the MTA report.
Rec_age0 <- c(107,142,75,61,112)*1e6

#mean(Rec_age0)
#sd(Rec_age0)

# fishing mortality at MSY (FMSY=F35%) - the target fishing mortality
# Table 2 - the updated 2021 MTA threshold fishing mortality reference point proxy.
F <- matrix(c(0.011, 0.048,	0.082, 0.079,	0.069, 0.066,	0.042, 0.015,
              0.007, 0.03,	0.063, 0.079,	0.073, 0.071,	0.043, 0.014,
              0.006, 0.028,	0.062, 0.086,	0.082, 0.08,	0.048, 0.016,
              0.008, 0.034,	0.082, 0.12,	0.116, 0.114,	0.069, 0.022,
              0.009, 0.039,	0.09,	 0.127,	0.122, 0.119,	0.072, 0.023),
            nrow = 8, ncol = 5, byrow = FALSE, dimnames = list(age = 0:7, year = 1:5))
F

# natural mortality
M <- 0.2

#total mortality
Z <- F + M
#mean(Z)
#sd(Z)
#log(sd(Z))

# Von Bertalanffy growth parameters for both male and female
Linf = 46.6
K = 0.15



## SIMULATE ABUNDANCE ####
pop <- sim_pop(iter = length(nsims), ages, years, Rec_age0, Z, Linf, K)

### EXTRACT N AT AGE MATRIX ####
Nage <- map_df(pop, ~pluck(., "N") |>
  #pop[[1]]$N |>
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

