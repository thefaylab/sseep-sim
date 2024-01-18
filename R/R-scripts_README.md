# sseep-sim / R

## Main

**01_make_grid.R** : create a grid of a given cell size that represents the NMFS bottom trawl survey footprint and areas for offshore wind energy development; to be used when generating distributions of simulated species abundances.

**02_simulate_populations.R** : simulate a given species-like population based on life history parameters and stock assessment estimates based on the SimSurvey package function `sim_abundance` and over a user specified number of iterations.

**03_append_distributions.R** : distribute the simulated species-like abundances across the bottom trawl survey grid and using predictions from a spatiotemporal generalized linear mixed model (GLMM) fit to the historical survey observations.

## R / Functions
**sim_pop_fn.R** : a function to loop through `sim_abundance` for a given number of iterations

**sim_stratmean_fn.R** : a function to calculate the random stratified mean abundance index from simulated survey data 

## R / Species Trials

Summer flounder

Atlantic mackerel

Scup

## R / Validation

**BTS_set_den.R**

**catch_rate_valid.R**
