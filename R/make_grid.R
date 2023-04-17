### created: 03/15/2023
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
library(sdmTMB)
library(marmap)
library(oce)
suppressPackageStartupMessages(library(tidyverse))
theme_set(theme_bw())
#source()

sdmtmb.dir <- "../sseep-analysis/sdmtmb"
sseep.dir <- "../sseep-analysis"
utm_proj <- "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

#### LOAD DATA ####
strata <- sf::st_read(here(sseep.dir, "gis", "NEFSC_BTS_AllStrata_Jun2022.shp")) %>%
  rename(STRATUM = "Strata_Num")

# read in summer flounder data
#sf_fall <- readRDS(here("sdmtmb", "data", "sumflounder_fall.rds"))
#sf_spring <-readRDS(here("sdmtmb", "data", "sumflounder_spring.rds"))

# read in best fit model
# m6_fall <- readRDS(here("sdmtmb", "model-outputs", "m6_fall.rds"))

#### DATA WRANGLING ####
# pull out unique values of strata for filtering
#sf_strat <- unique(sf_fall$STRATUM)
#sf_spr_strat <- unique(sf_spring$STRATUM)

# pull out depth values for adding later to our grid
# depths <- sf_fall |>
#   group_by(EST_YEAR, STRATUM) |>
#   summarise(AVGDEPTH = mean(AVGDEPTH)) #< mean()

# filter the full bts strata based on sf occurring strata
# fall_strat <- strata |>
#   filter(STRATUM %in% sf_strat)
#
# ggplot() +
#   geom_sf(data = fall_strat)

# x <- sf::st_coordinates(fall_strat) |> as.data.frame()
# y <- sdmTMB::add_utm_columns(x, c("X", "Y"), utm_names = c("a", "b"))

# IMPORTANT STEP: turn into UTMs first
# check my UTM choice
strata_utm <- sf::st_transform(strata, crs = 32618)

sf::st_crs(strata_utm)

ggplot() +
  geom_sf(data = strata_utm)

#### Make Fall Grid with SF ####
# choose a grid size in units of our polygon shape file
grid_spacing <- 10000 # 8000 m

# create grid over the bounding box of the polygon
full_grid <- sf::st_make_grid(
  strata_utm,
  cellsize = c(grid_spacing, grid_spacing),
  square = FALSE #creates hexagonal grid rather than square grid
) |>
  sf::st_sf()


# plot the grid
ggplot(full_grid) + geom_sf()

full_grid|>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  ggplot(aes(X, Y)) + geom_tile(width = grid_spacing, height = grid_spacing, colour = "grey40", fill = "white")

# subset our grid to cells that intersect our polygon:
intersected <- sf::st_intersects(full_grid, strata_utm)

selected_grid <- full_grid[lengths(intersected) > 0, ]

nrow(selected_grid) #4064

selected_grid|>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  ggplot(aes(X, Y)) + geom_tile(width = grid_spacing, height = grid_spacing, colour = "grey40", fill = "white")

# plot it
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = selected_grid, fill = NA)

join <- sf::st_join(selected_grid, strata_utm, largest = TRUE)
ggplot(join) + geom_sf()

# find how much of each grid cell is within the outer polygon:

# IMPORTANT STEP HERE!! take the union of your strata:
# strata_utm_union <- sf::st_union(strata_utm)
#
# ggplot(strata_utm_union) + geom_sf()
#
# overlap <- sf::st_intersection(selected_grid, strata_utm_union)
# nrow(overlap)
# nrow(selected_grid)
#
# ggplot(overlap) + geom_sf()

# calculate the area for index calculations later
# calculated_area <- sf::st_area(overlap)
# length(calculated_area)

# find the center points of each grid cell, extract the coordinates, and add the area values
coord <- join|>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  bind_cols(join) |>
  select(X, Y, STRATUM) #|>
  # st_as_sf()
  # mutate(across(c(X, Y), round, digits = 2)) |>
  #


ggplot(coord, aes(X, Y)) +
  geom_tile(width = grid_spacing, height = grid_spacing, colour = "grey10") +
  scale_fill_viridis_c() +
  coord_equal()
#
# ggplot(coord, aes(X, Y, fill = area)) +
#   geom_raster() +
#   scale_fill_viridis_c() +
#   coord_equal()

saveRDS(coord, here("data", "strata_coords.rds"))

#### BATHY ####

bts <- getNOAA.bathy(lon1 = -80, lon2 = -60,
                     lat1 = 32, lat2 = 46,
                     resolution = 0.25,
                     keep = TRUE)

grid_crs <- utm2lonlat(coord$X, coord$Y, zone = 18, hemisphere = "N")

#sf::st_crs(fall_strat_utm)

depths <- get.depth(bts, grid_crs, locator = FALSE)

grid <- bind_cols(coord, depths) |>
  select(X, Y, STRATUM, lon, lat, depth) |>
  mutate(cell = seq(length(STRATUM)),
         depths = abs(depth))
ggplot(grid, aes(X,Y))

#grid <- st_set_geometry(x = c(grid$lon, grid$lat), value = grid)

###Converting a dataframe to a raster in R
# coordinates(grid) <- ~ X + Y
# proj4string(grid) <- utm_proj
# grid <- spTransform(grid, crs(utm_proj))
# #Tell R that this is gridded:
# gridded(grid) = TRUE
#
#dat <- data.frame(cell = seq.int(length(grid)))
# survey_grid <- SpatialPixelsDataFrame(grid, dat)
# grid_ras <- raster(grid)

# plot(grid)
#
#
#
#
# as(grid, "SpatialPixels")

saveRDS(grid, here("data", "survey_grid.rds"))


