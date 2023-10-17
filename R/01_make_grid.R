### created: 03/15/2023
### last updated: 10/03/2023

# 01 - MAKE SURVEY GRID ####


## OBJECTIVE ####
# create a grid of the survey area
# identify a grid cell based on its occurrence within a given stratum
# add a wind ID code
# find the depth within each grid cell


### LOAD PACKAGES ####
library(stringr)
library(sf)
library(here)
library(marmap)
library(raster)
suppressPackageStartupMessages(library(tidyverse))
theme_set(theme_bw())
#source()

#sdmtmb.dir <- "../sseep-analysis/sdmtmb"
sseep.analysis <- "C:/Users/amiller7/Documents/cinar-osse/sseep-analysis"
utm_proj <- "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

### LOAD DATA ####
strata <- readRDS(here(sseep.analysis, "data", "rds", "active_strata.rds"))

# read in wind areas where they are one large polygon
wind_areas <- readRDS(here(sseep.analysis, "data", "rds", "wind_areas_062022", "merged_wind_areas_Jun2022.rds"))

# read in east coastline
east_coast <- sf::st_read(here("gis", "eastern_coast_UTM.shp"))

### DATA WRANGLE ####
# IMPORTANT STEP: turn into UTMs first
strata_utm <- sf::st_transform(strata, crs = 32618)
wind_utm <- sf::st_transform(wind_areas, crs = 32618)

# check that the crs transformation was agreeable
sf::st_crs(strata_utm)
sf::st_crs(wind_utm)

ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = wind_utm, fill = "lightblue", alpha = 0.5 )

# merge all strata for cell area calculation later in script
#strata_union <- sf::st_union(strata_utm)

## Make Grid ####
# choose a grid size in units of our polygon shape file
grid_spacing <- 10000 # 8000 m

# create grid over the bounding box of the polygon
full_grid <- sf::st_make_grid(
  strata_utm,
  cellsize = grid_spacing, #c(grid_spacing, grid_spacing),
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
# selected_grid <- sf::st_intersection(full_grid, strata_union)
# names(selected_grid)
#
# selected_grid <- selected_grid |>
#   rename(geometry = `sf..st_make_grid.strata_utm..cellsize...grid_spacing..square...FALSE.`) |>
#   mutate(calc_area = st_area(overlap_area),
#          cell_area_km = as.numeric(calc_area/1000000))


selected_grid <- full_grid[lengths(intersected) > 0, ]
#st_area(selected_grid)

#nrow(selected_grid) #2994

selected_grid |>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  ggplot(aes(X, Y)) + geom_tile(width = grid_spacing, height = grid_spacing, colour = "grey40", fill = "white")

# plot it
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = selected_grid, fill = NA)

### ADD STRATA ####
# find the intersection of each cell based on the largest overlap within a given stratum and attach that stratum data to the cell
join <- sf::st_join(selected_grid, strata_utm, largest = TRUE) |>
  mutate(cell = seq(1:length(geometry))) # add cell value which is referenced in SimSurvey functions

#overlap_area_join <- sf::st_join(overlap_area, strata_utm, largest = TRUE) |>
#  mutate(cell = seq(1:length(geometry)))

# plot
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = join, fill = NA)

### CREATE WIND ID ####
# combine multipolygons for faster intersect calculation
wind_utm_union <- sf::st_union(wind_utm)

# find the intersection of each cell based on its intersection within a given wind polygon; if intersects, value > 0
join$AREA_CODE <- sf::st_intersects(join, wind_utm_union) |>
  as.integer() |>
  dplyr::coalesce(2L) # find NAs within column and replace with the number 2 to represent the "outside wind area" ID

ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = join, aes(fill = AREA_CODE))

## EXTRACT BATHYMETRIC DATA ####
# use the marmap package to pull the bathymetric data from NOAA's database based on our survey area
bts <- getNOAA.bathy(lon1 = -80, lon2 = -60,
                     lat1 = 32, lat2 = 46,
                     resolution = 0.25,
                     keep = TRUE)

# convert the bathy object to an sf point object to manipulate with the other grid objects
bathy_grid <- marmap::as.xyz(bts) |>
  rename(lon = V1, lat = V2, depth = V3) |>
  sf::st_as_sf(coords = c("lon", "lat"))

# set the crs of the grid as the same as the survey strata because of the lat/long values
sf::st_crs(bathy_grid) <- sf::st_crs(strata)

# transform the lat/long values to utm values
bathy_grid <- sf::st_transform(bathy_grid, crs = 32618)

# spatial join the grid with the depth points to find all the depth points in a given cell
bathy_intersect <- sf::st_join(join, bathy_grid, left = FALSE) |> # inner spatial join
  group_by(cell) |>
  summarise(depth = mean(depth)) #|> # find the average depth in each cell
  #mutate(AREA_SqM = sf::st_area(cell))

# plot it
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = bathy_intersect, aes(fill = depth))


### REMOVE LAND CELLS IN THE GRID ####
# identify cells that are completely within the east coastline; if intersects, values > 0
land_cells <- sf::st_intersects(bathy_intersect, east_coast)

# remove the land cells
bathy_intersect <- bathy_intersect[lengths(land_cells) == 0, ]

# plot it
ggplot() +
  geom_sf(data = east_coast) +
  #geom_sf(data = strata_utm) +
  geom_sf(data = bathy_intersect, fill = NA)

## EXTRACT CENTROID COORDINATES OF GRID ####
# find the center points of each grid cell, extract the coordinates, and bind back cells and depths
grid_coords <- bathy_intersect |>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  bind_cols(bathy_intersect) |>
  mutate(across(c(X, Y), round, digits = 2)) |>
  dplyr::select(-geometry)

# plot it
ggplot(grid_coords, aes(X, Y, fill = depth)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_viridis_c() +
  coord_equal()

# ggplot() +
#   geom_sf(data = strata_utm, fill = "yellow") +
#   geom_tile(data = grid_coords, aes(X, Y), width = grid_spacing, height = grid_spacing, colour = NA, alpha = 0.5) +
#   coord_sf()

## FINAL MERGE ####
# create dataframe of cells and keep important data columns
join_df <- join |>
  sf::st_set_geometry(NULL) |>
  dplyr::select(STRATUM, cell, AREA_CODE)#, cell_area_km)

### A SIMPLE FEATURE ####
survey_grid_sf <- bathy_intersect |>
  left_join(join_df, by = "cell") |>
  mutate(AREA = case_when(          #add AREA_CODE for metadata
    AREA_CODE == 1 ~ "WIND",
    AREA_CODE == 2 ~ "OUTSIDE"),
    depth = depth * (-1),
    geometry = geometry/1000) |> # convert geometries to kms to match sdmtmb predictions and simsurvey units
  st_as_sf()

ggplot(survey_grid_sf) +
  geom_sf(aes(fill = as.factor(AREA), color = as.factor(AREA)))


### A DATAFRAME ####
# join merge data frames to combine all final columns,
survey_grid <- grid_coords |>
  left_join(join_df, by = "cell") |>
  mutate(AREA = case_when(          #add AREA_CODE for metadata
    AREA_CODE == 1 ~ "WIND",
    AREA_CODE == 2 ~ "OUTSIDE"),
    depth = depth * (-1),
    X = X/1000,
    Y = Y/1000) |>#,
    #cell_area_km = round(cell_area_km, 0)) |>    #multiply depths by -1 to match depth values used in BTS data
  as.data.frame()

# plot it
ggplot(survey_grid, aes(X, Y, fill = as.factor(AREA))) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_discrete() +
  coord_equal()

### SAVE GRID ####
saveRDS(survey_grid, here("data", "rds", "survey_grid_062022.rds"))

st_write(survey_grid_sf, here("gis", "survey_grid_062022.shp"), append = FALSE)

## CONVERT SF GRID FOR USE IN SIMSURVEY FNS ####
# to stars object
#test <- st_as_sf(survey_grid, coords = c("X", "Y")) |> st_rasterize() |> st_as_stars()
survey_grid_stars <- survey_grid_sf |>
  rename(strat = STRATUM) |>
  st_rasterize() |>
  st_as_stars()

saveRDS(survey_grid_stars, here("data", "survey_grid_stars_062022.rds"))
