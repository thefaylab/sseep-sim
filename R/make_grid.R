### created: 03/15/2023
### last updated: 04/24/2023

#### 01 - MAKE SURVEY GRID ####


# OBJECTIVE ####
# create a grid of the survey area
# identify a grid cell based on its occurrence within a given stratum
# add a wind ID code
# find the depth within each grid cell


## LOAD PACKAGES ####
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

#sdmtmb.dir <- "../sseep-analysis/sdmtmb"
#sseep.dir <- "../sseep-analysis"
utm_proj <- "+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## LOAD DATA ####
strata <- sf::st_read(here("gis", "NEFSC_BTS_AllStrata_Jun2022.shp")) %>%
  rename(STRATUM = "Strata_Num")

# read in wind areas where they are one large polygon
wind_areas <- sf::st_read(here("gis", "wind_areas_merge2023.shp"))

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
  geom_sf(data = wind_utm)

## Make Grid ####
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

# find the intersection of each cell based on the largest overlap within a given stratum and join that stratum data
join <- sf::st_join(selected_grid, strata_utm, largest = TRUE) |>
  mutate(cell = seq(1:length(geometry))) # add cell value which is referenced in SimSurvey functions

# plot
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = join, fill = NA)

## CREATE WIND ID ####
# combine multipolygons for faster intersect calculation
wind_utm_union <- sf::st_union(wind_utm)

# find the intersection of each cell based on its intersection within a given wind polygon; if intersects, value > 0
join$AREA_CODE <- st_intersects(join, wind_utm_union) |>
  as.integer() |>
  dplyr::coalesce(2L) # find NAs within column and replace with the number 2 to represent the "outside wind area" ID

### EXTRACT BATHYMETRIC DATA ####
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
st_crs(bathy_grid) <- st_crs(strata)

# transform the lat/long values to utm values
bathy_grid <- st_transform(bathy_grid, crs = 32618)

# spatial join the grid with the depth points to find all the depth points in a given cell
bathy_intersect <- sf::st_join(join, bathy_grid, left = FALSE) |> # inner spatial join
  group_by(cell) |>
  summarise(AVGDEPTH = mean(depth)) # find the average depth in each cell

# plot it
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = bathy_intersect, aes(fill = AVGDEPTH))


### REMOVE LAND CELLS IN THE GRID ####
# identify cells that are completely within the east coastline; if intersects, values > 0
land_cells <- st_within(bathy_intersect, east_coast)

# remove the land cells
bathy_intersect <- bathy_intersect[lengths(land_cells) == 0, ]

# plot it
ggplot() +
  geom_sf(data = east_coast) +
  #geom_sf(data = strata_utm) +
  geom_sf(data = bathy_intersect, fill = NA)

### EXTRACT CENTROID COORDINATES OF GRID ####
# find the center points of each grid cell, extract the coordinates, and bind back cells and depths
grid_coords <- bathy_intersect |>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  bind_cols(bathy_intersect) |>
  mutate(across(c(X, Y), round, digits = 2)) |>
  dplyr::select(-geometry)

# plot it
ggplot(grid_coords, aes(X, Y, fill = AVGDEPTH)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_viridis_c() +
  coord_equal()

ggplot() +
  geom_sf(data = strata_utm, fill = "yellow") +
  geom_tile(data = grid_coords, aes(X, Y), width = grid_spacing, height = grid_spacing, colour = NA, alpha = 0.5) +
  coord_sf()

### FINAL MERGE ####
# create dataframe of cells and keep important data columns
join_df <- join |>
  st_set_geometry(NULL) |>
  dplyr::select(STRATUM, cell, AREA_CODE)

# join merge data frames to combine all final columns,
survey_grid <- grid_coords |>
  left_join(join_df, by = "cell") |>
  mutate(AREA = case_when(          #add AREA_CODE for metadata
    AREA_CODE == 1 ~ "WIND",
    AREA_CODE == 2 ~ "OUTSIDE"),
    AVGDEPTH = AVGDEPTH * (-1))     #multiply depths by -1 to match depth values used in BTS data

# plot it
ggplot(survey_grid, aes(X, Y, fill = as.factor(AREA))) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_discrete() +
  coord_equal()

### SAVE GRID ####
saveRDS(survey_grid, here("data", "survey_grid.rds"))

## CONVERT GRID FOR USE IN SIMSURVEY FNS ####
# first to raster
#grid_ras <- raster::rasterFromXYZ(survey_grid, crs = utm_proj)

# then to RasterStack
#sdmtmb_stck <- stack(survey_grid)

# save the RasterStack file
#writeRaster(sdmtmb_ras, filename = here("data", "survey_grid.grd"))

# ggplot() +
#   geom_sf(data = strata_utm, fill = "yellow") +
#   geom_tile(data = coord, aes(X, Y), width = grid_spacing, height = grid_spacing, colour = NA, alpha = 0.5) +
#   geom_tile(data = sdmtmb_grid_fil, aes(X, Y), width = grid_spacing, height = grid_spacing, fill = "red") +
#   scale_fill_viridis_c() +
#   xlim(359969, 500000) +
#   ylim(3850000, 3950000) +
#   coord_sf()
