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
set.seed(380)
theme_set(theme_bw())

sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"

### DATA ####
# historical data
obs_catch <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data", "scup_fall.rds")) |>
  filter(EST_YEAR %in% c(2013:2014))


# the survey grid as dataframe
grid_xy <- readRDS(here("R", "scup", "data","survey_grid_Jun2022.rds")) |>
  rename(strat = STRATUM,
         x = X,
         y = Y) |>#, # rename to match column calls within SimSurvey functions
  #depth = AVGDEPTH) |>
  mutate(division = 1,
         depth = AVGDEPTH, #, # add division information
  x = x/1000, # convert to km
  y = y/1000) |>
  dplyr::select(x, y, cell, division, strat, depth) |>
  data.table::as.data.table()

# the survey Grid as stars object
grid_stars <- readRDS(here("R","scup","data",  "survey_grid_stars_062022.rds"))

# load spring model predictions
fall_preds <- readRDS(here(sseep.analysis, "sdmtmb", "scup", "data","fall_grid_preds.rds"))

# Recruitment at age 0 for the most recent year, 2013 and 2014; numbers were reported in thousands in the 66th SAW report.
Rec_age0 <- c(61, 112)*1e6

# fishing mortality at age for the most recent year, 2016 and 2017
ages <- as.character(0:7)
years <- as.character(1:2)
F <- matrix(c(0.008,0.034,0.082,0.120,0.116,0.114,0.069,0.022,
              0.009,0.039,0.090,0.127,0.122,0.119,0.072,0.023),
            nrow = 8, ncol = 2, byrow = FALSE, dimnames = list(age = 0:7, year = 1:2))
F
# natural mortality
M <- 0.2

#total mortality
Z <- F + M

# Von Bertalanffy growth parameters for both male and female
Linf = 46.6
K = 0.15


## SIMULATE ABUNDANCE ####
pop <- sim_abundance(ages = 0:7, years = 1:2,
                     R = sim_R(log_mean = log(Rec_age0), log_sd = 0.7, plot = TRUE),
                     Z = sim_Z(log_mean = log(Z), log_sd = 0.044, plot = TRUE),
                     N0 = sim_N0(N0 = "exp", plot = TRUE),
                     growth = sim_vonB(Linf = Linf, K = K, plot = TRUE))

pop$N0 <- c(60549,61129,91605,53072,35898,50528,46351,123923)*1000 # initial numbers at age from 2013
names(pop$N0) <- ages
pop$N <- array(c(60549,61129,91605,53072,35898,50528,46351,123923,
                 112436,49179,48375,69104,38540,26161,36895,134653)*1000,
               dim = c(8,2), dimnames = list(age = 0:7, year = 1:2))
#N_at_length <- convert_N(N_at_age = pop$N, lak = pop$sim_length(age = 0:7, length_age_key = TRUE))
#pop$N_at_length <- N_at_length

Nage <- tibble(age = 0:7,Nage = pop$N0)
plot_surface(pop, mat = "N")



## APPEND DISTRIBUTION ####

preds_2yr <- fall_preds |>
  filter(EST_YEAR %in% c(2013,2014)) |>
  rename(year = EST_YEAR, #rename to match column calls in SimSurvey
         strat = STRATUM) |>
  mutate(N_dist = exp(est),
         year = case_when( # change the year values based on their sequence
           year == 2013 ~ 1,
           year == 2014 ~ 2
          )) |>#,
  #cell = seq(1:length(N_dist))) |> # add cell # value
  dplyr::select(X,Y, year, N_dist, cell, strat, AVGDEPTH) |>
  data.table::as.data.table()



#plot it
ggplot(preds_2yr) +
  geom_tile(aes(X, Y, fill = N_dist), width = 10, height = 10) +
  scale_fill_viridis_c() +
  facet_wrap(~year) +
  labs(x = "Longitude", y = "Latitude", fill = "Number ") +
  theme_bw() +
  theme(legend.position = "bottom")

# predicted numbers at age based on simulated N0 and sdmtmb predictions
dist <- sdmTMB::replicate_df(preds_2yr, "age", c(0:7)) |> # replicate the predictions over each age, so there are distributions for each age
  left_join(Nage, by = "age") |> # join populations at age to the predictions
  group_by(year,age) |>
  rename(x = X,
         y = Y,
         depth = AVGDEPTH) |>
  mutate(P_i = N_dist/sum(N_dist), # probability of distribution
         N = Nage * P_i, # multiply probability by simulated numbers of age
         age = as.double(age),
         cell = as.double(cell),
         division = 1) |>
  dplyr::select(age, year, cell, N, x, y, strat, depth)

# forced numbers at age outside sdmTMB prediction area
no_dist <- dplyr::anti_join(grid_xy, preds_2yr, by = "cell") |>
  sdmTMB::replicate_df("age", c(0:7)) |>
  sdmTMB::replicate_df("year", c(1:2)) |>
  mutate(N = 0, #
         age = as.double(age),
         cell = as.double(cell),
         division = 1)

# bind distributions
full_dist <- bind_rows(no_dist, dist) |>
  as.data.table()


full_dist |>
  ggplot() +
  geom_tile(aes(x,y, fill = N), width = 10, height = 10) +
  scale_fill_viridis_c() + facet_wrap(~age)

# append  to abundance and rename list item
pop[[10]] <- grid_stars # use stars object of  grid
names(pop)[[10]] <- "grid"

# append to the simulated abundance object
pop[[11]] <- grid_xy
names(pop)[[11]] <- "grid_xy"

# append to the simulated abundance object
pop[[12]] <- full_dist
names(pop)[[12]] <- "sp_N"


# plot numbers at age
ggplot(pop$sp_N) +
  geom_tile(aes(x, y, fill = N), width = 10, height = 10) +
  scale_fill_viridis_c(trans = "sqrt") +
  facet_wrap(~age, scales = "free_y") +
  labs(x = "Longitude", y = "Latitude", fill = "Number of Scup") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45))

#saveRDS(pop, here("R", "scup", "data", "forced-pop-dist.rds"))

#-----------------------------------------------------------------
#function to set sel to half
sim_logistic_custom <- function(k = 2, x0 = 3, plot = TRUE) {
    function(x = NULL) {
      y <- (1 / (1 + exp(-k * (x - x0))))*0.5
      if (plot) plot(x, y, type = "b")
      y
    }
  }
#-----------------------------------------------------------------
#
#
#
## SIMULATE SURVEY ####
surv_dat <- pop |>
  sim_survey(n_sims = 20,
             trawl_dim = c(2.7, 0.014),
             q = sim_logistic_custom(k = 2, x0 = 2.5), #new function modifies selectivity to catch half of the fish by age
             set_den = 0.001,
             min_sets = 3,
             #lengths_cap = # max # of lengths to collect per set
             #age_cap = 8, # max # of ages to sample per length group - all 8?
             age_sampling = "stratified",
             #age_length_group = , # length group bin size (cm) -
             age_space_group = "set", # spatial scale of stratified age sampling - ages taken each tow?
             resample_cells = TRUE)

#saveRDS(pop, here("data", "rds", "catch-rate-valid", "sim100-survey.rds"))

# plot it
data <- surv_dat$setdet
wind_areas <- readRDS(here("R", "scup","data", "all_wind_areas_Jun2022.rds"))
strata <- readRDS(here("R", "scup","data", "active_strata.rds"))

wind_areas_utm <- st_transform(wind_areas, crs = 32618)
strata_utm <- st_transform(strata, crs = 32618)

#simulated catch rate
ggplot() +
  geom_sf(data = strata_utm, fill = NA) +
  geom_sf(data = wind_areas_utm, fill = "lightblue") +
  coord_sf() +
  geom_point(data = data |> filter(year == 1), aes(x*1000,y*1000, size = n, color = n)) +
  scale_color_continuous(type = "viridis") +
  scale_size_continuous(guide = "none") +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", color = "Number") +
  theme(legend.position = "bottom", panel.spacing =  unit(1, "lines"), axis.text.x = element_text(angle = 45, vjust = 0.5))



h1<-ggplot() +
  geom_histogram(data = obs_catch |> filter(EST_YEAR == 2013), aes(EXPCATCHNUM), bins = 30) +
  xlim(NA,1000) #+ ylim(NA,4)

h2<- ggplot() +
  geom_histogram(data = data |> filter(year == 1), aes(n), bins = 30) +
  xlim(NA,1000) #+ ylim(NA,100)

grid.arrange(h1, h2, ncol=1)

## VALIDATE CATCH RATE ####
obs_strat <- unique(obs_catch$STRATUM)

obs_resamps <- obs_catch |>
  rep_slice_sample(reps = 20, replace = TRUE, prop = 1)

#saveRDS(pop, here("data", "rds", "catch-rate-valid", "resamp100-obs-data.rds"))


# obs_catch |> group_by(EST_YEAR) |> mutate(TOWID = str_c(STRATUM, CRUISE6, STATION)) |> summarise(towct = length(unique(TOWID)), mu_catch = (sum(EXPCATCHNUM)/towct), sd = sd(EXPCATCHNUM)) |> mutate(sim = 1) |> rename(year = EST_YEAR)
#
# obs_catch |> mutate(TOWID = str_c(STRATUM, CRUISE6, STATION)) |> summarise(towct = length(unique(TOWID)), mu_catch = (sum(EXPCATCHNUM)/towct), sd = sd(EXPCATCHNUM))
#
# data |> filter(strat %in% obs_strat) |> group_by(year, sim) |> summarise(towct = length(unique(set)), mu_catch = (sum(n)/towct), sd = sd(n))

# data wrangle
obs_data <- obs_resamps |>
  mutate(TOWID = str_c(STRATUM, CRUISE6, STATION),
         TYPE = "Observed",
         EST_YEAR = case_when(
           EST_YEAR == 2013 ~ 1,
           EST_YEAR == 2014 ~ 2
         )) |>
  select(replicate, TOWID, EXPCATCHNUM, EST_YEAR, TYPE, AVGDEPTH) |>
  rename(n = EXPCATCHNUM,
         year = EST_YEAR,
         depth = AVGDEPTH)

sim_data <- data |>
  dplyr::filter(strat %in% obs_strat) |>
  mutate(TOWID = as.character(seq(set)),
         TYPE = "Simulated") |>
  select(sim, TOWID, n, year, TYPE, depth) |>
  rename(replicate = sim)

test_data <- bind_rows(obs_data, sim_data)

### MEAN TEST ####
# mu_test <- t.test(n ~ TYPE,  data = test_data, alternative = "two.sided", mu = 0, conf.level = 0.95, paired = FALSE)
#
# mu_test_df <- map_dfr(mu_test, ~pluck(., 1)) |>
#   bind_cols(map(mu_test, ~pluck(., 2))) |>
#   janitor::clean_names() |>
#   rename(lower.ci = conf_int_4,
#          upper.ci = conf_int_11,
#          simulated_mu = estimate_12,
#          observed_mu = estimate_5,
#          formula = data_name,
#          tstat = statistic,
#          df = parameter) |>
#   select(!c(alternative, null_value))
# mu_test_df <- mu_test_df[, c(8,7, 5, 10, 6, 4, 9, 1:3)]
#
# kable(mu_test_df |> select(!method), caption = mu_test_df$method) |> kable_styling()

### PROPORTION OF ZEROS ####
# proportions of zeros in sim is very diff from proportion of zeros of real data - why is that?
prop_zero <- test_data |>
  group_by(TYPE, replicate, year) |>
  nest() |>
  mutate(zero = map(data, ~filter(., n == 0) |> nrow()),
         total = map(data, ~nrow(.))) |>
  select(!data) |>
  unnest(cols = c(zero, total)) |>
  summarise(prop = round((zero/total)*100, 0))

saveRDS(prop_zero, here("data", "rds", "catch-rate-valid", "prop-zeros.rds"))

# kable(prop_zero, caption = "Proportion of zero catch rates") |> kable_styling(full_width = FALSE)

ggplot(prop_zero) +
  geom_boxplot(aes(x = as.factor(year), y = prop, color = TYPE)) +
  labs(x = "Year", y = "Proportions of zeros (%)", subtitle = "Distribution of proportion of zeros") + ylim(0,90)
# distributions don't look that different but the mean is larger in the simulated data

ggsave("dist-prop-zeros.png", device = "png", path = here("outputs", "scup", "catch-rate-valid"), width = 8, height = 5)

### DISTRIBUTION PLOTS ####
# catch rates
ggplot(test_data) +
  geom_histogram(aes(y = n, fill = TYPE)) +
  facet_wrap(year~TYPE) +
  labs(y = "Catch per tow", x = "Frequency of catch rate occurrence", subtitle = "Distribution of catch rates")
# distributions don't look that different but the mean is larger in the simulated data

ggsave("histogram-catch-rates.png", device = "png", path = here("outputs", "scup", "catch-rate-valid"), width = 8, height = 5)

# ggplot(test_data) +
#   geom_boxplot(aes(x = as.factor(year), y = n, color= TYPE)) +
#   labs(x = "Year", y = "Frequency of catch rate occurrence", subtitle = "Distribution of catch rates")
# # distributions don't look that different but the mean is larger in the simulated data
#
# ggsave("boxplot-catch-rates.png", device = "png", path = here("outputs", "scup", "catch-rate-valid"), width = 8, height = 5)

write.csv(test_data, file = "test.data.csv")
# average catch rate
test_data2 <- test_data |>
 group_by(replicate,TYPE, year) |>
  summarise(mu = mean(n))

(test_data_plot <- test_data2 |>
  ggplot() +
  geom_boxplot(aes(x = as.factor(year), y = mu, color= TYPE)) + ylim(0,NA) +
  labs(x = "Year", y = "Mean catch", subtitle = ""))
# distributions don't look that different but the mean is larger in the simulated data

ggsave("boxplot-mu_catch-rates.png", device = "png", path = here("outputs", "scup", "catch-rate-valid"), width = 8, height = 5)


# tows with respect to depth
# ggplot(test_data) +
#   geom_histogram(aes(x = depth, fill = TYPE)) +
#   facet_wrap(~TYPE)

ggplot(test_data) +
  geom_point(aes(x = depth, y = n)) +
  facet_wrap(~TYPE) +
  labs(x = "Depth (m)", y = "Catch rate (N/tow)", subtitle = "Distribution of catch rates with respect to depth")

ggsave("depth-catch-rates.png", device = "png", path = here("outputs", "scup", "catch-rate-valid"), width = 8, height = 5)


