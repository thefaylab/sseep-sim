### created: 04/05/2023
### last updated:

####  ####

###################
#### OBJECTIVE ####
###################
#


####################


#### LOAD PACKAGES ####
library(stringr)
library(sf)
library(patchwork)
library(here)
library(raster)
library(SimSurvey)
library(marmap)
library(oce)
suppressPackageStartupMessages(library(tidyverse))
theme_set(theme_bw())
#source()

sdmtmb.dir <- "../sseep-analysis/sdmtmb"

#### LOAD DATA ####
# load summer flounder data
sumflounder <- readRDS(file = here(sdmtmb.dir, "data", "sumflounder.rds"))

# load spring model
m7_spring <- readRDS(file = here(sdmtmb.dir, "model-outputs", "m7_spring.rds"))

# load spring model predictions
spring_preds1 <- readRDS(file = here(sdmtmb.dir, "data", "spring_preds.rds"))

sdmtmb_grid <- readRDS(file = here(sdmtmb.dir, "data", "survey_grid.rds")) |>
  rename(strat = STRATUM) |>
  mutate(division = 1)


# filter distribution predictions for year 2021
preds2019 <- spring_preds1 |>
  filter(EST_YEAR == 2019) |>
  rename(N_dist = est,
         year = EST_YEAR) |>
  dplyr::select(c(X,Y, year, N_dist)) |>
  mutate(year = 1,
         cell = seq(1:length(N_dist))) |>
  data.table::as.data.table()

#### STOCK ASSESSMENT DATA ####
## 66th SAW Report
Rec_age0 <- c(48689, 60598, 44582, 33088, 28416) #most recent 5 years
mean(Rec_age0)
sd(Rec_age0)

F <- c(0.340, 0.286, 0.331, 0.414, 0.419) #most recent 5 years
M <- 0.2

Z <- F + M

mean(Z)
sd(Z)
log(sd(Z))

#### SIMULATE ABUNDANCE ####
set.seed(8675309)
pop <- sim_abundance(ages = 0:7, years = 1,
                     R = sim_R(log_mean = log(mean(Rec_age0)), log_sd = 0.001, plot = TRUE),
                     Z = sim_Z(log_mean = log(mean(Z)), log_sd = 0.001, plot = TRUE),
                     growth = sim_vonB(Linf = 83.6, K = 0.14, plot = TRUE)) #66th SAW Report ASFMC
Nage <- tibble(age = 0:7, pop$N0)
plot_surface(pop, mat = "N")


#### CONVERT GRID ####
sdmtmb_ras <- rasterFromXYZ(sdmtmb_grid, crs = 32618) # able to convert with raster pkg 3.6-11, errors occuring with newer verions
writeRaster(sdmtmb_ras, filename = here("data", "survey_grid.grd"))

# append  to abundance and rename list item
pop <- append(pop,sdmtmb_ras)
names(pop)[[10]] <- "grid"

##### EXPAND GRID TO INCLUDE AGE AND YEAR ####
grid_xy <- sdmtmb_grid |>
  dplyr::select(X, Y, depth, division, strat) |>
  rename(x = X, y = Y) |>
  mutate(cell = seq(1:length(depth)))
pop <- append(pop, list(grid_xy))
names(pop)[[11]] <- "grid_xy"

##### APPEND sdmTMB DISTRIBUTION ####
# calculate numbers at age in each cell based on the year distribution

dist <- sdmTMB::replicate_df(preds2019, "age", c(0:7)) |>
  left_join(Nage, by = "age") |>
  rename(Nage = "pop$N0") |>
  mutate(N_dist = exp(N_dist),
         P_i = N_dist/sum(N_dist),
         N = Nage * P_i) |>
  dplyr::select(age, year, cell, N, X, Y)
pop <- append(pop, list(dist))
names(pop)[[12]] <- "sp_N"

ggplot(pop$sp_N) + geom_tile(aes(X, Y, fill = N), width = 10, height = 10) + scale_fill_viridis_c() + facet_wrap(~age)

# str(pop_dist$grid_xy)
# str(pop_dist$sp_N)
# nrow(preds2019)
# nrow(pop_dist$grid_xy)

#### SIMULATE DISTRIBUTION TRIALS ####
# survey_grid <- readRDS(file = here("data", "survey_grid.rds")) |>
#   mutate(X = X/1000,
#          Y = Y/1000)
# survey_grid <- rasterFromXYZ(survey_grid, crs = 32618)
#
# #with temporary grid
# g <- make_grid(res = c(10, 10), n_div = 4, strat_splits = 2,
#                shelf_depth = 100, shelf_width = 80, depth_range = c(10, 375))
#
# pop_dist1 <- sim_distribution(pop,grid = g,
#                               ays_covar = sim_ays_covar(range = 89.5),
#                               depth_par = sim_parabola(mu = 71.79,
#                                                        sigma = 52.92))
# # with sdmtmb grid
# pop_dist2 <- sim_distribution(pop, grid = survey_grid,
#                    ays_covar = sim_ays_covar(range = 89.5),
#                    depth_par = sim_parabola(mu = 71.79,
#                                             sigma = 52.92))
#
# plot_trend(pop_dist, sum_ages = pop_dist$ages, col = viridis::viridis(1))
# plot_distribution(pop_dist1, ages = 0:7, years = 1, type = "heatmap",scale = "natural")
# plot_distribution(pop_dist2, ages = 0:7, years = 1, type = "heatmap",scale = "natural")
# plot_distribution(pop_dist, ages = 0:7, years = 1, type = "contour")

