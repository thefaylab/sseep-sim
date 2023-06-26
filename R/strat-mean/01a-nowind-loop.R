### created: 04/27/2022
### last updated: 06/26/2023

# 01a - STRATIFIED CALCULATIONS: WITH WIND AREAS INCLUDED ####


## OBJECTIVE ####
# calculate stratified means for each species, strata, and year combination based on the historical time series data. The analysis contained herein calculates the stratified means for the full time series without any change to represent the base/control scenario. No tows/observations have been subtracted from this dataset



## LOAD PACKAGES ####
library(stringr)
library(sf)
library(patchwork)
library(here)
suppressPackageStartupMessages(library(tidyverse))

#sdmtmb.dir <- "../sseep-analysis/sdmtmb"
sseep.analysis <- "../sseep-analysis"

## LOAD DATA ####
#source(here(sseep.analysis, "R", "StratMeanFXs_v2.R"))

# dataset created from `03b-spatial-filter-data.R` here(sseep_analysis,"tidy-data"). Contains complete observations for each species and unique tow filtered based on 95% cumulative distribution of biomass.
data <- readRDS(here("data", "rds", "95filtered_complete_bts.rds"))

# species dataframe for adding to final dataset
species <- readRDS(here(sseep.analysis, "data", "rds", "95filtered-species.rds"))

#
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "strata_wts.rds"))
BTSArea <- sum(strata_wts$Area_SqNm)


# CALCULATE INDIVIDUAL MEANS AND VARIANCES ####
# calculate individual means and variances for each combinations of species, year, and stratum
nowind_means <- data |>
  group_by(STRATUM, EST_YEAR, SVSPP, SEASON) |> #, GEO_AREA) %>%
  summarise(towct = length(unique(STATION)), # calculate unique tows
            mu = sum(EXPCATCHNUM)/towct, # find the average biomass based on unique tows rather than observations to avoid potential duplication
            var = ifelse(towct == 1, 0, # if the tow count equals 1, then variance about the mean should be 0
                         sum((EXPCATCHNUM - mu)^2)/(towct - 1))) %>% # if tow count does not equal 1, then find the variance of biomass
  left_join(strata_wts, by = "STRATUM") %>% # add each stratum area and relative weight to the dataset based on STRATUM number
  mutate(wt_mu = Area_SqNm * mu, # part one of the stratified mean formula
         wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) # part one of the stratified variance formula


### save the data
saveRDS(nowind_means, here("data", "rds", "indiv-num-mu_included.rds"))



# COMPLETE STRATIFIED MEAN AND VARIANCE CALCULATIONS ####
# calculate stratified means and variances for each combinations of species and year based on individual stratum means and variances
nowind_stratmu <- nowind_means %>%
  group_by(SVSPP, EST_YEAR, SEASON) %>% #, GEO_AREA) %>%
  summarise(stratmu = (sum(wt_mu)) / BTSArea, # part two of the stratified mean formula
            stratvar = sum(wt_var)) %>% # part two of the stratified variance formula
  mutate(TYPE = paste("With Wind Included")) # paste identifying information of means and variances for joining and plotting in later scripts


### save the data
saveRDS(nowind_stratmu, here("data", "rds", "num-strat-mu_included.rds"))


# tow count by area
# data %>%
#   mutate(code = str_c(STRATUM, CRUISE6, STATION)) %>%
#   #group_by(SVSPP, SEASON, AREA) %>%
#   summarise(towct = length(code)) #%>%
#
# data %>%
#   mutate(code = str_c(STRATUM, CRUISE6, STATION)) %>%
#   #group_by(SVSPP, SEASON, AREA) %>%
#   summarise(towct = length(unique(code))) #%>%
#
# data %>%
#   mutate(code = str_c(STRATUM, CRUISE6, STATION)) %>%
#   group_by(SEASON, AREA) %>%
#   summarise(towct = length(unique(code))) #%>%
#   #left_join(specieslookup, by = "SVSPP")
#
#
##################################################

#### CALCULATE INDIVIDUAL MEANS AND VARIANCES BY YEAR####
# calculate individual means and variances for each combinations of species, year, and stratum
nowind_means_yr <- data %>%
  #mutate(tow_code = str_c(STRATUM, CRUISE6, STATION)) %>%
  #left_join(geounits, by = "code") %>%
  #rename(STRATUM = STRATUM.x,
  #       STATION = STATION.x) %>%
  group_by(STRATUM, EST_YEAR, SVSPP, COMNAME) %>% #, GEO_AREA) %>%
  summarise(towct = length(unique(STATION)), # calculate unique tows
            mu = sum(EXPCATCHNUM)/towct, # find the average biomass based on unique tows rather than observations to avoid potential duplication
            var = ifelse(towct == 1, 0, # if the tow count equals 1, then variance about the mean should be 0
                         sum((EXPCATCHNUM - mu)^2)/(towct - 1))) %>% # if tow count does not equal 1, then find the variance of biomass
  left_join(strata_wts, by = "STRATUM") %>% # add each stratum area and relative weight to the dataset based on STRATUM number
  mutate(wt_mu = Area_SqNm * mu, # part one of the stratified mean formula
         wt_var = ((((RelWt)^2) * var) / towct) * (1 - (towct / Area_SqNm))) # part one of the stratified variance formula


### save the data
saveRDS(nowind_means_yr, here("data", "rds", "indiv-mu_included.rds"))



#### COMPLETE STRATIFIED MEAN AND VARIANCE CALCULATIONS BY YEAR####
# calculate stratified means and variances for each combinations of species and year based on individual stratum means and variances
nowind_stratmu_yr <- nowind_means_yr %>%
  group_by(SVSPP,COMNAME, EST_YEAR) %>% #, GEO_AREA) %>%
  summarise(stratmu = (sum(wt_mu)) / BTSArea, # part two of the stratified mean formula
            stratvar = sum(wt_var)) %>% # part two of the stratified variance formula
  mutate(TYPE = paste("With Wind Included")) # paste identifying information of means and variances for joining and plotting in later scripts


### save the data
saveRDS(nowind_stratmu_yr, here("data", "rds", "strat-mu_included.rds"))


##
(strat_yr <- nowind_stratmu_yr |>
    filter(SVSPP %in% c(141, 105, 106, 121,15,131,32,72,503,22,23,24,25,26,27,28, 103, 143)))
write_csv(strat_yr, file = "strat_yr.csv")




