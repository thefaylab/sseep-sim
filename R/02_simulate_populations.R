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
source(here("R", "sim_pop_fn.R"))
# random seed number for generating different populations but for reproducibility of simulation
set.seed(131)


### DATA SET UP ####
# location of where to save data
pop.dat <- here("data", "rds", "pops")
Nage.dat <- here("data", "rds", "Nages")

### name of species to be simulated
species <- "sumflounder"

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

# Recruitment at age 0 for the most recent 5 years, 2015-2019; numbers were reported in thousands in the MTA report.
Rec_age0 <- c(28416, 33088, 44582, 60598, 48689)*1000

#mean(Rec_age0)
#sd(Rec_age0)

# fishing mortality at MSY (FMSY=F35%) - the target fishing mortality
# Table 2 - the updated 2021 MTA threshold fishing mortality reference point proxy.
F <- 0.422

# natural mortality
M <- 0.25

#total mortality
Z <- F + M
#mean(Z)
#sd(Z)
#log(sd(Z))

# Von Bertalanffy growth parameters for both male and female
Linf = 83.6
K = 0.14



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


## SAVE THE DATA ####
saveRDS(pop, here(pop.dat, str_c(species, length(nsims), "pop.rds", sep = "_")))
saveRDS(Nage, here(Nage.dat, str_c(species, length(nsims), "Nage.rds", sep = "_")))






