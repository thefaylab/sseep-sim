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
spring_preds <- readRDS(file = here(sdmtmb.dir, "data", "spring_preds.rds"))

# filter distribution predictions for year 2021
preds2019 <- spring_preds |>
  filter(EST_YEAR == 2019) |>
  rename(N = est,
         year = EST_YEAR) |>
  select(c(X,Y, year, N)) |>
  mutate(year = 1) |>
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
round(pop$N)
plot_surface(pop, mat = "N")


#### CONVERT GRID ####
sdmtmb_ras <- rasterFromXYZ(sdmtmb_grid, crs = 32618)

# append  to abundance and rename list item
pop <- append(pop,sdmtmb_ras)
names(pop)[[10]] <- "grid"

#### EXPAND GRID TO INCLUDE AGE AND YEAR ####
sdmtmb_egrid <- sdmtmb_grid
i <- rep(seq(nrow(sdmtmb_egrid)), times = length(0:7)) # replicate over number of ages
a <- rep(0:7, each = nrow(sdmtmb_egrid))
sdmtmb_egrid <- sdmtmb_egrid[i, ]
sdmtmb_egrid$age <- a
sdmtmb_egrid <- sdmtmb_egrid[,!3:5]
# i <- rep(seq(nrow(sdmtmb_egrid)), times = 1)
# y <- rep(1, each = nrow(sdmtmb_egrid))
# sdmtmb_egrid <- sdmtmb_egrid[i, ]
# sdmtmb_egrid$year <- y
sdmtmb_egrid$cell <- seq(1:nrow(sdmtmb_egrid))

# grid_xy$strat
# xy <- sdmtmb_egrid[,c("X", "Y")]

#error <- ays_covar(x = xy, ages = seq(0:7), years = 1, cells = sdmtmb_egrid$cell)

pop <- append(pop, list(sdmtmb_egrid))
names(pop)[[11]] <- "grid_xy"
pop <- append(pop, list(preds2019))
names(pop)[[12]] <- "sp_N"

# str(pop_dist$grid_xy)
# str(pop_dist$sp_N)
# nrow(preds2019)
# nrow(pop_dist$grid_xy)

#### SIMULATE DISTRIBUTION ####
# temporary grid
g <- make_grid(res = c(10, 10), n_div = 4, strat_splits = 1,
               shelf_depth = 100, shelf_width = 80, depth_range = c(10, 375))

pop_dist <- sim_abundance(ages = 0:7, years = 1,
                          R = sim_R(log_mean = log(mean(Rec_age0)), log_sd = 0.001, plot = TRUE),
                          Z = sim_Z(log_mean = log(mean(Z)), log_sd = 0.001, plot = TRUE),
                          growth = sim_vonB(Linf = 83.6, K = 0.14, plot = TRUE)) |>
  sim_distribution(grid = g,
                   ays_covar = sim_ays_covar(range = 89.5),
                   depth_par = sim_parabola(mu = 71.79,
                                            sigma = 52.92))

plot_trend(pop_dist, sum_ages = pop_dist$ages, col = viridis::viridis(1))
plot_distribution(pop_dist, ages = 0:7, years = 1, type = "heatmap",scale = "natural")
plot_distribution(pop_dist, ages = 0:7, years = 1, type = "contour")

