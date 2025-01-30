### created:
### updated:
# VALIDATE CATCH RATES FOR Scup ####


### LOAD PACKAGES ####
library(SimSurvey)
suppressPackageStartupMessages(library(tidyverse))
library(data.table)
library(patchwork)
library(here)
library(infer)
library(sf)
library(gridExtra)
library(kableExtra)
set.seed(380)
theme_set(theme_bw())

sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"

### DATA ####
# historical data
obs_catch_fall <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data", "scup_fall.rds")) |>
  filter(EST_YEAR %in% c(2010,2011,2012,2013,2014))

obs_catch_fall_full <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data", "scup_fall.rds")) |>
  filter(EST_YEAR %in% c(2009:2016,2018,2019,2021))
obs_catch_fall_full$ID <- seq_len(nrow(obs_catch_fall_full))

obs_catch_fall |>
  group_by(EST_YEAR) |>
  summarise(mu_num = mean(EXPCATCHNUM),
            mu_wt = mean(EXPCATCHWT))


obs_catch_fall_full |> group_by(EST_YEAR) |> summarise(tows = length(unique(TOWID)))

ggplot(obs_catch_fall_full, aes(x = ID, y = AVGDEPTH)) +
  geom_point(color = "blue") + # Add points
  labs(title = "Depth per Tow - Data",
       x = "Tow ID",
       y = "Depth (m)") +
  theme_minimal()


ggplot(obs_catch_fall) +
  geom_point(aes(x = AVGDEPTH, y = EXPCATCHWT)) +
  facet_wrap(~YEAR, ncol = 5) +
  labs(x = "Depth (m)", y = "Catch (weight) per tow", subtitle = "")


# the survey grid as dataframe
grid_xy_m2 <- readRDS(here("data", "rds", "survey_grid_all_122024.rds")) |>
  rename(strat = STRATUM,
         x = X,
         y = Y,
         depth = median_2) |>#, # rename to match column calls within SimSurvey functions
  mutate(division = 1) |>#, # add division information
  dplyr::select(x,y,cell,depth, strat, AREA_CODE,division) |>
  data.table::as.data.table() |>
  drop_na()

# as stars object
grid_stars_2 <- readRDS(here("data",  "survey_grid_all_stars_122024.rds")) |>
  rename(depth = mean_2) |>
  dplyr::select(cell,depth, strat, AREA_CODE)



summary(grid_xy_1$depth)
summary(grid_xy_2$depth)
summary(grid_xy_m1$depth)
summary(grid_xy_m2$depth)

# the survey Grid as stars object
#grid_stars <- readRDS(here("R","scup","data",  "survey_grid_stars_062022.rds"))
grid_stars_2

# load spring model predictions
fall_preds_1 <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data","fall_grid_preds_tw_1.rds"))
fall_preds_2 <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data","fall_grid_preds_tw_2.rds"))
fall_preds_m1 <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data","fall_grid_preds_tw_m1.rds"))
fall_preds_m2 <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data","fall_grid_preds_tw_m2.rds"))


fall_preds_2 <- fall_preds_2 |>
  dplyr::select(!c(median_1, mean_1, n1, Cell_Area, median_2,mean_2, n2, New_area, Perc_reduction, Perc_Reduction))

source(file = here(sseep.analysis,"sdmtmb", "R", "plot_fns.R"))
plot_preds(fall_preds_m1, exp(est)) + scale_fill_viridis_c(trans = "sqrt", option = "H") +
  labs(fill = "Biomass estimates") +
  theme(legend.position = "right")

# Recruitment at age 0 for the most recent year, 2013 and 2014; numbers were reported in thousands in the 66th SAW report.
Rec_age0 <- c(107,142,75,61,112)*1e6

# fishing mortality at age for the most recent year
ages <- as.character(0:7)
years <- as.character(1:5)
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

# Von Bertalanffy growth parameters for both male and female
Linf = 46.6
K = 0.15


## SIMULATE ABUNDANCE ####
pop <- sim_abundance(ages = 0:7, years = 1:5,
                     R = sim_R(log_mean = log(Rec_age0), log_sd = 0.001, plot = TRUE),
                     Z = sim_Z(log_mean = log(Z), log_sd = 0.001, plot = TRUE),
                     N0 = sim_N0(N0 = "exp", plot = TRUE),
                     growth = sim_vonB(Linf = Linf, K = K, plot = TRUE))

pop$N0 <- c(107141,79663,177355,106495,85639,59115,29172,78993)*1000 # initial numbers at age from 2010
names(pop$N0) <- ages

pop$N <- array(c(107141, 79663,	 177355, 106495, 85639,	59115, 29172,	78993,
                 141523, 86802,	 62159,	 88502,	 80586,	65439, 45302,	86619,
                 75149,	 115086, 68981,	 47781,	 66981,	61340, 49896,	105457,
                 60549,	 61129,	 91605,	 53072,	 35898,	50528, 46351,	123923,
                 112436, 49179,	 48375,	 69104,	 38540,	26161, 36895,	134653) *1000,
               dim = c(8,5), dimnames = list(age = 0:7, year = 1:5))

Nage <- as.data.frame.table(pop$N) |>
  arrange(age) |>
  rename(Nage = Freq) |>
  mutate(sim = 1) |>
  relocate(sim) |>
  as_tibble() |>
  mutate(sim = as.integer(sim),
         year = as.integer(year),
         age = as.numeric(levels(age))[age])|>
  mutate(age = as.integer(age))

plot_surface(pop, mat = "N")


## APPEND DISTRIBUTION ####

preds_5yr_m1 <- fall_preds_m1 |>
  filter(EST_YEAR %in% c(2010,2011,2012,2013,2014)) |>
  rename(year = EST_YEAR, #rename to match column calls in SimSurvey
         strat = STRATUM) |>
  mutate(N_dist = exp(est),
         year = case_when( # change the year values based on their sequence
           year == 2010 ~ 1,
           year == 2011 ~ 2,
           year == 2012 ~ 3,
           year == 2013 ~ 4,
           year == 2014 ~ 5
          )) |>#,
  #cell = seq(1:length(N_dist))) |> # add cell # value
  dplyr::select(X,Y, year, N_dist, cell, strat, AVGDEPTH) |>
  data.table::as.data.table()


#plot it
ggplot(preds_5yr_m1) +
  geom_tile(aes(X, Y, fill = N_dist), width = 10, height = 10) +
  scale_fill_viridis_c() +
  facet_wrap(~year) +
  labs(x = "Longitude", y = "Latitude", fill = "Biomass") +
  theme_bw() +
  theme(legend.position = "bottom")

summary(preds_5yr_1$AVGDEPTH)

# predicted numbers at age based on simulated N0 and sdmtmb predictions
dist_m2 <- sdmTMB::replicate_df(preds_5yr_m2, "age", c(0:7)) |> # replicate the predictions over each age, so there are distributions for each age
  left_join(Nage, by = c("age","year")) |> # join populations at age and year to the predictions
  rename(x = X,
         y = Y,
         depth = AVGDEPTH) |>
  group_by(year, age) |>
  mutate(P_i = N_dist/sum(N_dist), # probability of distribution
         N = Nage * P_i, # multiply probability by simulated numbers of age
         age = as.double(age),
         cell = as.double(cell),
         division = 1) |>
  dplyr::select(age, year, cell, N, x, y, strat, depth)

summary(dist_2$depth)

# forced numbers at age outside sdmTMB prediction area
no_dist_1 <- dplyr::anti_join(grid_xy_1, preds_5yr_1, by = "cell") |>
  sdmTMB::replicate_df("age", c(0:7)) |>
  sdmTMB::replicate_df("year", c(1:5)) |>
  mutate(N = 0, #
         age = as.double(age),
         cell = as.double(cell),
         division = 1)

summary(no_dist_1$depth)

# bind distributions
full_dist_1 <- bind_rows(no_dist_1, dist_1) |>
  as.data.table()


full_dist_2 |>
  ggplot() +
  geom_tile(aes(x,y, fill = N), width = 10, height = 10) +
  scale_fill_viridis_c() + facet_wrap(~age)


pop_2 <- pop
# append  to abundance and rename list item
pop_2[[10]] <- grid_stars_2 # use stars object of  grid
names(pop_2)[[10]] <- "grid"

# append to the simulated abundance object
pop_2[[11]] <- grid_xy_2
names(pop_2)[[11]] <- "grid_xy"

# append to the simulated abundance object
pop_2[[12]] <- full_dist_2
names(pop_2)[[12]] <- "sp_N"


# plot numbers at age
ggplot(pop_m2$sp_N) +
  geom_tile(aes(x, y, fill = N), width = 10, height = 10) +
  scale_fill_viridis_c(trans = "sqrt") +
  facet_wrap(~age, scales = "free_y") +
  labs(x = "Longitude", y = "Latitude", fill = "Biomass of Scup") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45))

#saveRDS(pop_m2, here("R", "scup", "data", "forced-pop-dist_1.rds"))
pop_1 <- readRDS(here("R", "scup", "data", "forced-pop-dist_1.rds"))
pop_2 <- readRDS(here("R", "scup", "data", "forced-pop-dist_2.rds"))
pop_m1 <- readRDS(here("R", "scup", "data", "forced-pop-dist_m1.rds"))
pop_m2 <- readRDS(here("R", "scup", "data", "forced-pop-dist_m2.rds"))

summary(pop_1$sp_N$depth)
summary(pop_2$sp_N$depth)

summary(pop_m1$sp_N$depth)
summary(pop_m2$sp_N$depth)


#-----------------------------------------------------------------
#function to set selectivity

source(here("R/selectivity_fns.R"))

#-----------------------------------------------------------------
## SIMULATE SURVEY ####
surv_dat_m2 <- pop_m2 |>
  sim_survey(n_sims = 25,
             trawl_dim = c(2.7, 0.014),
             q = force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1),
             set_den = 0.0001,
             min_sets = 2,
             #lengths_cap = # max # of lengths to collect per set
             #age_cap = 8, # max # of ages to sample per length group - all 8?
             age_sampling = "stratified",
             #age_length_group = , # length group bin size (cm) -
             age_space_group = "set", # spatial scale of stratified age sampling - ages taken each tow?
             resample_cells = TRUE)



plot_survey(surv_dat_m1, which_year = 1)
plot_survey(surv_dat_m2, which_year = 1)
# plot_survey(surv_dat, which_year = 2)
# plot_survey(surv_dat, which_year = 3)
# plot_survey(surv_dat, which_year = 4)
# plot_survey(surv_dat, which_year = 5)



# plot it
data_m2 <- surv_dat_m2$setdet

ggplot(data_m2, aes(x = sim, y = depth)) +
  geom_point(color = "orange") + # Add points
  labs(title = "Depth per Tow - Sim Survey",
       x = "Tow ID",
       y = "Depth (m)") +
  theme_minimal()

wind_areas <- readRDS(here("R", "scup","data", "all_wind_areas_Jun2022.rds"))
strata <- readRDS(here("R", "scup","data", "active_strata.rds"))

wind_areas_utm <- st_transform(wind_areas, crs = 32618)
strata_utm <- st_transform(strata, crs = 32618)

#simulated catch rate
ggplot() +
  geom_sf(data = strata_utm, fill = NA) +
  geom_sf(data = wind_areas_utm, fill = "lightblue") +
  coord_sf() +
  geom_point(data = data_m2, aes(x*1000,y*1000, size = n, color = n)) +
  facet_wrap(~year) +
  scale_color_continuous(type = "viridis") +
  scale_size_continuous(guide = "none") +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", color = "Number", title="Sim Catch Rate", subtitle="median filtered") +
  theme(legend.position = "right", panel.spacing =  unit(1, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5))



## VALIDATE CATCH RATE ####
obs_catch_fall_area <- obs_catch_fall |>
  left_join(strata, by="STRATUM")

obs_strat <- unique(obs_catch_fall_area$STRATUM)


# obs_resamps_old <- obs_catch_fall_area |>
# rep_slice_sample(reps = 100, replace = TRUE, prop = 1) #resample with proportionality to the strata
#sampling spatially proportionally to each strata to preserve the space
#saveRDS(obs_catch_fall_area,"obs_catch_fall_area.rds")
#

#--------------------------------------------
# Calculate the total area
total_area <- sum(obs_catch_fall_area$Area_SqNm)

# Calculate the probability for each row based on the stratum's area
obs_catch_fall_area <- obs_catch_fall_area %>%
  group_by(STRATUM) %>%
  mutate(prob = Area_SqNm / total_area) %>%
  ungroup()

# Perform the proportional resampling
obs_resamps_new <- obs_catch_fall_area %>%
  rep_slice_sample(reps = 100, replace = TRUE, weight_by = prob, prop = 1)

#----------------------------------------------------------

# data wrangle
obs_data <- obs_resamps_new |>
  mutate(TOWID = str_c(STRATUM, CRUISE6, STATION),
         TYPE = "Observed",
         EST_YEAR = case_when(
           EST_YEAR == 2010 ~ 1,
           EST_YEAR == 2011 ~ 2,
           EST_YEAR == 2012 ~ 3,
           EST_YEAR == 2013 ~ 4,
           EST_YEAR == 2014 ~ 5
         )) |>
   dplyr::select(replicate, TOWID, EXPCATCHNUM, EST_YEAR, TYPE, AVGDEPTH) |>
  rename(n = EXPCATCHNUM,
         year = EST_YEAR,
         depth = AVGDEPTH)

obs_catch <- obs_catch_fall |>
  mutate(TOWID = str_c(STRATUM, CRUISE6, STATION),
         TYPE = "Observed",
         replicate = NA,
         EST_YEAR = case_when(
           EST_YEAR == 2010 ~ 1,
           EST_YEAR == 2011 ~ 2,
           EST_YEAR == 2012 ~ 3,
           EST_YEAR == 2013 ~ 4,
           EST_YEAR == 2014 ~ 5
         )) |>
  dplyr::select(replicate, TOWID, EXPCATCHNUM, EST_YEAR, TYPE, AVGDEPTH) |>
  rename(n = EXPCATCHNUM,
         year = EST_YEAR,
         depth = AVGDEPTH)




###----SURVEY FIXED LOCATIONS-------------------------------------------
#load a previous survey object
### DATA SET UP ####
# data locations
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
survdat <- here("data", "rds", "survdat")

species <- "scup"
season <- "fall"
nsims <- 1:2
nsurveys <- 25
ages <- 0:7
years <- 1:5

survdat_sq_1 <- readRDS(here(survdat, str_c(species, season, length(nsims), "sims", nsurveys, "sq-surv-dat_1.rds", sep = "_")))

survdat_sq_1 <- survdat_sq_1[[1]]

fixed_grid_1 <- survdat_sq_1$grid_xy

# Access the grid from the pop object
pop_grid_1 <- pop_1$grid_xy

# Subset the pop grid to include only the cells from the fixed grid
# cell is the column that contains the grid cell IDs in both pop and fixed_locations
grid_fixed_1 <- pop_grid_1[pop_grid_1$cell %in% fixed_grid_1$cell, ]

# Update pop object to use only the fixed grid cells
pop_fixed_1 <- pop_1
pop_fixed_1$grid_xy <- grid_fixed_1

# survey with fixed locations will need a different minimum set parameter, without replacement
surv_dat_fixed_m1 <- pop_fixed_m1 |>
  sim_survey(n_sims = 100,
             trawl_dim = c(2.7, 0.014),
             q = force_sim_logistic(k = -0.66, x0 = -1.14, plot = TRUE, force_age = TRUE, age = 0, force_sel = 1),
             set_den = 0.0001,
             min_sets = 3,
             age_sampling = "stratified",
             age_space_group = "set",
             resample_cells = FALSE)



saveRDS(surv_dat_fixed_m1, here("R", "scup", "data", "surv_dat_fixed_m1.rds"))
surv_dat_fixed_1<- readRDS(here("R", "scup","data", "surv_dat_fixed_m1.rds"))
surv_dat_fixed_2<- readRDS(here("R", "scup","data", "surv_dat_fixed_m2.rds"))
surv_dat_fixed_m1<- readRDS(here("R", "scup","data", "surv_dat_fixed_m1.rds"))
surv_dat_fixed_m2<- readRDS(here("R", "scup","data", "surv_dat_fixed_m2.rds"))



# Extract sampled locations from surv_dat_fixed
sampled_locations_1 <- surv_dat_fixed_1$grid$cell  # Modify based on the structure of surv_dat_fixed


#Data frame for fixed and sampled locations
fixed_df_1 <- data.frame(x = fixed_grid_1$x, y = fixed_grid_1$y)
sampled_df_1 <- data.frame(x = surv_dat_fixed_1$grid_xy$x, y = surv_dat_fixed_1$grid_xy$y)

# Plot the fixed and sampled locations
ggplot() +
  geom_point(data = fixed_df_1, aes(x = x, y = y), color = "blue", size = 3, alpha = 0.5) +
  geom_point(data = sampled_df_1, aes(x = x, y = y), color = "red", size = 2, alpha = 0.7) +
  labs(title = "Fixed vs Sampled Locations",
       x = "Long", y = "Lat")


surv_fixed_setdet_1 <- surv_dat_fixed_1$setdet
surv_fixed_setdet_2 <- surv_dat_fixed_2$setdet
surv_fixed_setdet_m1 <- surv_dat_fixed_m1$setdet
surv_fixed_setdet_m2 <- surv_dat_fixed_m2$setdet

summary(surv_fixed_setdet_1$depth)
summary(surv_fixed_setdet_2$depth)
summary(surv_fixed_setdet_m1$depth)
summary(surv_fixed_setdet_m2$depth)

# plot it

#data2 <- surv_dat_fixed$setdet

#simulated catch rate
ggplot() +
  geom_sf(data = strata_utm, fill = NA) +
  geom_sf(data = wind_areas_utm, fill = "lightblue") +
  coord_sf() +
  geom_point(data = data2, aes(x*1000,y*1000, size = n, color = n)) +
  facet_wrap(~year) +
  scale_color_continuous(type = "viridis") +
  scale_size_continuous(guide = "none") +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", color = "Number") +
  theme(legend.position = "right", panel.spacing =  unit(1, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5))



fordepth <- data2 |>
  dplyr::filter(strat %in% obs_strat) |>
  mutate(TOWID = as.character(seq(set)),
         TYPE = "Simulated") |>
  dplyr::select(sim, TOWID, n, year, TYPE, depth, strat, cell, x, y) |>
  rename(replicate = sim)

fordepth |> filter(depth >100) |>
ggplot() +
  geom_sf(data = strata_utm, fill = NA) +
  #geom_sf(data = wind_areas_utm, fill = "lightblue") +
  coord_sf() +
  geom_point(aes(x*1000,y*1000, size = n, color = n)) +
 # facet_wrap(~year) +
  scale_color_continuous(type = "viridis") +
  scale_size_continuous(guide = "none") +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", color = "Number") +
  theme(legend.position = "right", panel.spacing =  unit(1, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5))

depth_out<-fordepth #|> filter(depth >100)

depth_out <- depth_out %>%
  arrange(cell) %>%
  dplyr::select(cell, depth, strat) |>
  rename(depth_mean=depth) |>
  distinct()

unique(depth_out$cell)
table(depth_out$cell,depth_out$strat)




###---------------------------------------------------------------------
# data wrangle with survey w/fixed locations
# obs data is the same


sim_data2 <- data2 |>
  dplyr::filter(strat %in% obs_strat) |>
  mutate(TOWID = as.character(seq(set)),
         TYPE = "Simulated") |>
  dplyr::select(sim, TOWID, n, year, TYPE, depth) |>
  rename(replicate = sim)

test_data2 <- bind_rows(obs_catch, sim_data2)

### PROPORTION OF ZEROS ####
# proportions of zeros in sim is very diff from proportion of zeros of real data - why is that?
prop_zero <- test_data |>
  group_by(TYPE, replicate, year) |>
  nest() |>
  mutate(zero = map(data, ~filter(., n == 0) |> nrow()),
         total = map(data, ~nrow(.))) |>
  dplyr::select(!data) |>
  unnest(cols = c(zero, total)) |>
  summarise(prop = round((zero/total)*100, 0))

prop_zero2 <- test_data2 |>
  group_by(TYPE, replicate, year) |>
  nest() |>
  mutate(zero = map(data, ~filter(., n == 0) |> nrow()),
         total = map(data, ~nrow(.))) |>
  dplyr::select(!data) |>
  unnest(cols = c(zero, total)) |>
  summarise(prop = round((zero/total)*100, 0))


ggplot(prop_zero) +
  geom_boxplot(aes(x = as.factor(year), y = prop, color = TYPE)) +
 labs(x = "Year", y = "Proportions of zeros (%)", subtitle = "Distribution of proportion of zeros \n Random Location Survey")

ggplot(prop_zero2) +
  geom_boxplot(aes(x = as.factor(year), y = prop, color = TYPE)) +
  labs(x = "Year", y = "Proportions of zeros (%)", subtitle = "Distribution of proportion of zeros \n Fixed location Survey")


### DISTRIBUTION PLOTS ####
# catch rates
# ggplot(test_data) +
#   geom_histogram(aes(y = n, fill = TYPE)) +
#   facet_wrap(year~TYPE) +
#   labs(y = "Catch per tow", x = "Frequency of catch rate occurrence", subtitle = "Distribution of catch rates")


ggplot(test_data) +
  geom_boxplot(aes(x = as.factor(year), y = n, color= TYPE)) +
  labs(x = "Year", y = "Frequency of catch rate occurrence", subtitle = "Distribution of catch rates \n Random Location Survey")

ggplot(test_data2) +
  geom_boxplot(aes(x = as.factor(year), y = n, color= TYPE)) +
  labs(x = "Year", y = "Frequency of catch rate occurrence", subtitle = "Distribution of catch rates \n Fixed location Survey")

#write.csv(test_data, file = "test.data.csv")
# average catch rate
catch_rate <- test_data |>
  group_by(replicate,TYPE, year) |>
  summarise(mu = mean(n))

(catch_rate_plot <- catch_rate |>
    ggplot() +
    geom_boxplot(aes(x = as.factor(year), y = mu, color= TYPE)) + ylim(0,NA) +
    labs(x = "Year", y = "Mean catch", subtitle = "Catch rate \n Random Location Survey"))

catch_rate2 <- test_data2 |>
  group_by(replicate,TYPE, year) |>
  summarise(mu = mean(n))

(catch_rate_plot2 <- catch_rate2 |>
    ggplot() +
    geom_boxplot(aes(x = as.factor(year), y = mu, color= TYPE)) + ylim(0,NA) +
    labs(x = "Year", y = "Mean catch", subtitle = "Catch rate \n Fixed Location Survey"))


test_data |>
  group_by(TYPE, year) |>
  summarise(mu = mean(n))

test_data2 |>
  group_by(TYPE, year) |>
  summarise(mu = mean(n))


data |>
  group_by(year, sim) |>
  summarise(tows = length(unique(set))) |>
  summarise(mean(tows))

data2 |>
  group_by(year, sim) |>
  summarise(tows = length(unique(set))) |>
  summarise(mean(tows))


obs_catch_fall |> group_by(EST_YEAR) |> summarise(tows = length(unique(TOWID)))

sum<-obs_catch_fall_area |> group_by(EST_YEAR) |> summarise(tows = length(unique(TOWID)))
sum1<-obs_catch_fall_area |> group_by(EST_YEAR,STRATUM) |> summarise(tows = length(unique(TOWID)))


##Strat mean
source(here("R", "sim_stratmean_fn.R"))
source(here(sseep.analysis,"R","StratMeanFXs_v2.R"))
strata_wts <- readRDS(here(sseep.analysis,"data", "rds", "active_strata_wts.rds"))
strata_wts2 <-  strata_wts |> rename(strat = STRATUM)


(strat_mean_OC <- obs_resamps_new |>
    mutate(year = YEAR,
           TYPE = "Observed",
           year = case_when(
             year == 2010 ~ 1,
             year == 2011 ~ 2,
             year == 2012 ~ 3,
             year == 2013 ~ 4,
             year == 2014 ~ 5
           )) |>
    group_by(replicate, year,TYPE) |>
    nest() |>
    mutate(strat_mu = map(data, ~stratified.mean(.,strata_wts, "number"))) |>
    unnest(cols = strat_mu) |>
    arrange(year))




(strat_mean_SC <- data |>
    mutate(TYPE = "Simulated")|>
    group_by(sim, year, TYPE) |>
    nest() |>
    mutate(strat_mu = map(data, ~sim_stratmean2(.,strata_wts2))) |>
    unnest(cols = strat_mu) |>
    rename(replicate = sim))


strat_mean_data <- bind_rows(strat_mean_OC, strat_mean_SC)


(strat_mean_p <- strat_mean_data |>
    ggplot() +
    geom_boxplot(aes(x = as.factor(year), y = stratmu, color= TYPE)) + ylim(0,NA) +
    labs(x = "Year", y = "Strat Mean catch", subtitle = ""))



ggplot(test_data) +
  geom_point(aes(x = depth, y = n)) +
  facet_wrap(~TYPE) +
  labs(x = "Depth (m)", y = "Catch rate (N/tow)", subtitle = "Distribution of catch rates with respect to depth")


ggplot(test_data2) +
  geom_point(aes(x = depth, y = n)) +
  facet_wrap(~TYPE) +
  labs(x = "Depth (m)", y = "Catch rate (N/tow)", subtitle = "Distribution of catch rates with respect to depth - fixed loc")







