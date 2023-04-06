### created: 04/05/2023
### last updated: 

#### LOAD PACKAGES ####
library(stringr)
library(sf)
library(patchwork)
library(here)
library(raster)
library(sdmTMB)
library(marmap)
library(oce)
library(sf)
library(raster)
library(SimSurvey)
library(dplyr)
library(ggplot2)
library(plotly)
library(terra)
library(ncdf4)
suppressPackageStartupMessages(library(tidyverse))
theme_set(theme_bw())
#source()

#### LOAD DATA ####
strata <- sf::st_read(here("gis", "NEFSC_BTS_AllStrata_Jun2022.shp")) %>%
  rename(STRATUM = "Strata_Num")

plot(strata)
class(strata)

strata_poly <- st_transform(strata, "WGS84")
plot(strata_poly)

utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"
strata_poly_utm <- strata_poly %>% st_transform(utm_proj)



# extract bathymetric data for survey area
bathy_bts <- getNOAA.bathy(lon1 = -80, lon2 = -60,
                           lat1 = 32, lat2 = 46, 
                           resolution = 0.25, 
                           keep = TRUE)
bathy_bts <- rast(bathy_bts)

strat_bathy <- raster::mask(bathy_bts, strata_poly) # cells within strata

## Strata-by-strata depths
strat_no <- unique(strata_poly$STRATUM)
means <- numeric(length(strat_no))
mins <- numeric(length(strat_no))
maxs <- numeric(length(strat_no))
names(means) <- names(mins) <- names(maxs) <- strat_no

for (s in strat_no) {
  
  temp <- raster::mask(strat_bathy, strata_poly[strata_poly$STRATUM == s, ])
  
  cat(s, "\n")
  cat("mean:", mean(temp[], na.rm = TRUE), "\n")
  cat("min:", min(temp[], na.rm = TRUE), "\n")
  cat("max:", max(temp[], na.rm = TRUE), "\n")
  
  means[as.character(s)] <- mean(temp[], na.rm = TRUE)
  mins[as.character(s)] <- min(temp[], na.rm = TRUE)
  maxs[as.character(s)] <- max(temp[], na.rm = TRUE)
  
}

depth_by_strata <- data.frame(strat = strat_no,
                              mean = means,
                              min = mins, max = maxs)

## Survey area
survey_area <- sum(sf::st_area(strata_poly_utm))
survey_area
# 330978.1 [km^2]

## Area of the stratum
strat_sums <- strata_poly_utm %>%
  mutate(area = sf::st_area(strata_poly_utm)) %>%
  group_by(STRATUM) %>%
  summarize(strat_area = sum(area))
range(strat_sums$strat_area)
# Units: [km^2]
# 55.33499 13886.26696

## Number of strata
strata <- length(unique(strata_poly_utm$STRATUM))
strata
# 176

## Mean area
mean_area <- mean(strat_sums$strat_area)
mean_area
# 1880.557 [km^2]

## Range depth in the survey area
depth_range <- range(strat_bathy[], na.rm = TRUE)
depth_range
# -1727.978   156.731

## Mean depth in the survey area
depth_mean <- mean(strat_bathy[], na.rm = TRUE)
depth_mean
# -94.85613

depth_median <- median(strat_bathy[], na.rm = TRUE)
depth_median
# -67.93962

## Determines x_y_range
sqrt(survey_area)/2
# 287.6535 [km]

## TO DO: Modify make_grid arguments to create a similar survey area
## shelf depth/width modifies mean grid depth, start breaks/splits modifies
## strata (length), method determines slope of shelf
library(stars)
library(bezier)
grid <- make_grid(x_range = c(-204, 204),
                  y_range = c(-204, 204),
                  res = c(3.5, 3.5),
                  shelf_depth = 100,
                  shelf_width = 100,
                  depth_range = c(100, 1800),
                  n_div = 1,
                  strat_breaks = seq(0, 2000, by = 50),
                  strat_splits = 2,
                  method = "bezier")

#tibble::lst(survey_area, strata, mean_area, depth_range, depth_mean, depth_median)


#grid2 = st_xy2sfc(grid, as_points = TRUE)
#grid3 <- st_as_sf(grid, as_points = TRUE, merge = FALSE)

prod(res(grid)) * ncell(grid)
#
length(unique(grid$strat))
#
mean(table(values(grid$strat)) * prod(res(grid)))
#
range(values(grid$depth), na.rm = TRUE)
#
mean(values(grid$depth), na.rm = TRUE)
#
median(values(grid$depth), na.rm = TRUE)
#

xyz <- data.frame(rasterToPoints(grid$depth))
plot_ly(data = xyz, x = ~x, y = ~-depth) %>% add_lines()

plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")

grid_dat <- data.frame(raster::rasterToPoints(grid))

strat_depths <- grid_dat %>%
  group_by(strat) %>%
  summarize(mean = mean(depth), min = min(depth), max = max(depth))

## Compare real vs simulated depths
real_depths <- data.frame(depth_by_strata)
real_depths <- abs(real_depths)
sim_depths <- data.frame(strat_depths)
all_depths <- rbind(real_depths, sim_depths)

all_depths <- all_depths %>% mutate(grid = factor(ifelse(strat > 300, "real", "sim")))

all_depths %>%
  filter(!is.na(grid)) %>%
  ggplot(aes(x=mean, color=grid, fill=grid)) +
  geom_histogram(position="dodge", bins=10, alpha = 1) + theme_bw()

plot_grid(grid)

