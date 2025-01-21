##SANDBOX - DEPTH ANALYSIS

# 01 - MAKE SURVEY GRID WITH MEDIANS CALCULATED BASED ON ALL THE DATA AND FILTERED FOR LESS THAN 400M DEPTH####


## OBJECTIVE ####
# create a grid of the survey area
# identify a grid cell based on its occurrence within a given stratum
# add a wind ID code
# find the depth within each grid cell

# loading the workspace
#load("saveworkspace_grid.RData")

### LOAD PACKAGES ####
library(stringr)
library(sf)
library(here)
library(marmap)
library(raster)
library(stars)
suppressPackageStartupMessages(library(tidyverse))
theme_set(theme_bw())

sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
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

## Make Grid ####
# choose a grid size in units of our polygon shape file
grid_spacing <- 10000

# create grid over the bounding box of the polygon
full_grid <- sf::st_make_grid(
  strata_utm,
  cellsize = grid_spacing, #c(grid_spacing, grid_spacing),
 # what = "polygons",
  square = FALSE) |> #creates hexagonal grid rather than square grid
  sf::st_sf()

# Add an area column (in square meters by default)
full_grid <- full_grid |>
  mutate(Cell_Area = as.numeric(sf::st_area(sf..st_make_grid.strata_utm..cellsize...grid_spacing..square...FALSE.)))  # Convert area to numeric if needed and change name because the new version of the sf package create this longer column names


# plot the grid
ggplot(full_grid) + geom_sf()

full_grid|>
   sf::st_centroid() |>
   sf::st_coordinates() |>
   as_tibble(.name_repair = "minimal") |>
   ggplot() +
   geom_sf(data = full_grid, fill = "white", color = "grey40") #+ # Hexagonal grid



# subset our grid to cells that intersect our polygon:
intersected <- sf::st_intersects(full_grid, strata_utm)
selected_grid <- full_grid[lengths(intersected) > 0, ]

# plot it
ggplot() +
  geom_sf(data = selected_grid, fill = "white", color = "grey40") #+ # Hexagonal grid


ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = selected_grid, fill = NA)

### ADD STRATA ####
# find the intersection of each cell based on the largest overlap within a given stratum and attach that stratum data to the cell

join <- sf::st_join(selected_grid, strata_utm, largest = TRUE) |>
  mutate(cell = row_number()) # add cell value which is referenced in SimSurvey functions

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

# Perform the spatial join first
joined_bathy_data <- sf::st_join(join, bathy_grid, left = FALSE)

cells_check_mean <- c(337,338,817,379)
cells_check_median <- c(337,338,817,379,888,704,339)

filter(grid_coords, cell %in% cells_check_mean) |>
ggplot(aes(X, Y)) +
  geom_tile(data=grid_coords, width = grid_spacing, height = grid_spacing, fill = "grey80") +
  geom_tile(width = grid_spacing, height = grid_spacing, fill="black") +
#  scale_fill_viridis_c() +
  coord_equal()


filter(grid_coords, cell %in% cells_check_median) |>
  ggplot(aes(X, Y)) +
  geom_tile(data=grid_coords, width = grid_spacing, height = grid_spacing, fill = "grey80") +
  geom_tile(width = grid_spacing, height = grid_spacing, fill="black") +
  #  scale_fill_viridis_c() +
  coord_equal()


# Calculate the required medians and counts
result_all_depths <- joined_bathy_data |>
  relocate(Cell_Area, .after = cell) |>
  group_by(cell) |>  # Group by cell
  summarise(
    median_1 = median(depth, na.rm = TRUE),  # Median for all depths
    mean_1 = mean(depth, na.rm = TRUE), # Mean for all depths
    n1 = n(),                               # Count for all depths
    Cell_Area = first(Cell_Area),         # Include cell area
    geometry = st_geometry(first(sf..st_make_grid.strata_utm..cellsize...grid_spacing..square...FALSE.)), # Retain geometry (sf property)
    .groups = "drop") |>
  st_as_sf()  # Ensure it remains an sf object

# Calculate the medians and counts for filtered data <400m
result_filtered_depths <- joined_bathy_data |>
  group_by(cell) |>  # Group by cell
  summarise(
    median_2 = median(depth[depth >= -400], na.rm = TRUE), # Median for depths <400m
    mean_2 = mean(depth[depth >= -400], na.rm = TRUE), # Mean for depths <400m
    n2 = sum(depth >= -400, na.rm = TRUE),
    geometry = st_geometry(first(sf..st_make_grid.strata_utm..cellsize...grid_spacing..square...FALSE.)), # Retain geometry (sf property)
    .groups = "drop") |>
  st_as_sf()  # Ensure it remains an sf object


#join both data sets
final_result <- st_join(
  result_all_depths,
  result_filtered_depths |> dplyr::select(-geometry),  # Drop geometry from the second dataset
  join = st_equals                              # Use exact spatial matches
) |>
  mutate(
    New_area = ifelse(n1 == 0, NA, (n2 / n1) * Cell_Area),
    Perc_reduction = ((Cell_Area - New_area)/Cell_Area)*100
  )


#final_result <- sf::st_transform(final_result, crs = 32618)

#convert into data frame
final_result_df <- as.data.frame(final_result)

#save
write.csv(final_result_df, "final_result_with_area_proportion.csv", row.names = FALSE)

### REMOVE LAND CELLS IN THE GRID ####
land <- sf::st_intersects(final_result, east_coast)
depth_intersect <- final_result[lengths(land) == 0, ]


##REMORE CELL N2 = 0 (perc reduction 100%) and other arrangements
bathy_intersect <- depth_intersect |>
 # filter(n2 >= 1) |>
  rename(cell = cell.x) |>
  dplyr::select(-cell.y)


# plot it
ggplot() +
  geom_sf(data = strata_utm) +
  geom_sf(data = bathy_intersect, aes(fill = median_1))



## EXTRACT CENTROID COORDINATES OF GRID ####
# find the center points of each grid cell, extract the coordinates, and bind back cells and depths
grid_coords <- bathy_intersect |>
  sf::st_centroid() |>
  sf::st_coordinates() |>
  as_tibble() |>
  bind_cols(bathy_intersect) |>
  mutate(across(c(X, Y), \(x) round(x, digits = 2))) |>
  dplyr::select(-geometry)

# plot it
# ggplot(result_mean_depths, aes(X, Y, fill = median_1)) +
#   geom_tile(width = grid_spacing, height = grid_spacing) +
#   scale_fill_viridis_c() +
#   coord_equal()

ggplot(grid_coords, aes(X, Y, fill = median_1)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_viridis_c() +
  coord_equal()

ggplot(grid_coords, aes(X, Y, fill = median_2)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_viridis_c() +
  coord_equal()

ggplot(grid_coords, aes(X, Y, fill = mean_1)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_viridis_c() +
  coord_equal()

ggplot(grid_coords, aes(X, Y, fill = mean_2)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_viridis_c() +
  coord_equal()



grid_coords <- grid_coords |>
   mutate(
      Perc_Reduction = case_when(
        Perc_reduction > 95 ~ ">95",                          # Greater than 95
        Perc_reduction > 70 & Perc_reduction <= 95 ~ "70-95", # Between 70 and 95
        Perc_reduction > 30 & Perc_reduction <= 70 ~ "30-70", # Between 30 and 70
        Perc_reduction <= 30 ~ "<30",                         # Lower than or equal to 30
      TRUE ~ "Unknown"))                                      # Catch any edge cases



# Ensure Perc_Reduction is a factor with the desired order
grid_coords <- grid_coords |>
  mutate(
    Perc_Reduction = factor(
      Perc_Reduction,
      levels = c("<30", "30-70", "70-95", ">95")))  # Specify the order


# Plot with custom legend order
ggplot(grid_coords, aes(X, Y, fill = Perc_Reduction)) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_manual(values = c("<30" = "forestgreen","30-70" = "gold","70-95" = "orange2",">95" = "firebrick2")) +
  coord_equal() +
  labs(fill = "Reduction Category",
  title = "Grid Percent Reduction Categories") +
  theme_minimal()


write.csv(grid_coords, "grid_coords_df.csv", row.names = FALSE)


## FINAL MERGE ####
# create dataframe of cells and keep important data columns
join_df <- join |>
  sf::st_set_geometry(NULL) |>
  dplyr::select(STRATUM, cell, AREA_CODE)#, cell_area_km)

### A SIMPLE FEATURE ####
survey_grid_all_sf <- bathy_intersect |>
  left_join(join_df, by = "cell") |>
  mutate(AREA = case_when(          #add AREA_CODE for metadata
    AREA_CODE == 1 ~ "WIND",
    AREA_CODE == 2 ~ "OUTSIDE"),
    median_1 = median_1 * (-1),
    median_2 = median_2 * (-1),
    mean_1 = mean_1 * (-1),
    mean_2 = mean_2 * (-1),
   geometry = geometry/1000) |> # convert geometries to kms to match sdmtmb predictions and simsurvey units
  st_as_sf()

ggplot(survey_grid_all_sf) +
  geom_sf(aes(fill = as.factor(AREA), color = as.factor(AREA)))

save.image("saveworkspace_grid.RData")

### A DATAFRAME ####
# join merge data frames to combine all final columns,
survey_grid_all <- grid_coords |>
  left_join(join_df, by = "cell") |>
  mutate(AREA = case_when(
    AREA_CODE == 1 ~ "WIND", #add AREA_CODE for metadata
    AREA_CODE == 2 ~ "OUTSIDE"),
    median_1 = median_1 * (-1),
    median_2 = median_2 * (-1),
    mean_1 = mean_1 * (-1),
    mean_2 = mean_2 * (-1),
    X = X/1000,
    Y = Y/1000) |>#,
    #cell_area_km = round(cell_area_km, 0)) |>    #multiply depths by -1 to match depth values used in BTS data
  as.data.frame()

# plot it
ggplot(survey_grid_all, aes(X, Y, fill = as.factor(AREA))) +
  geom_tile(width = grid_spacing, height = grid_spacing) +
  scale_fill_discrete() +
  coord_equal()

### SAVE GRID ####
saveRDS(survey_grid_all, here("data", "rds", "survey_grid_all_122024.rds")) #adds median and mean columns
st_write(survey_grid_all_sf, here("gis", "survey_grid_all_122024.shp"), append = FALSE)

## CONVERT SF GRID FOR USE IN SIMSURVEY FNS ####
# to stars object
survey_grid_all_stars <- survey_grid_all_sf |>
  rename(strat = STRATUM) |>
  st_rasterize() |>
  st_as_stars()

saveRDS(survey_grid_all_stars, here("data", "survey_grid_all_stars_122024.rds"))



###############################################################################################
###OBSERVED LOCATIONS###

# Define analysis directory and load data
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"

complete <- readRDS(here(sseep.analysis, "data", "rds", "completed_bts_data.rds")) %>%
  filter(EST_YEAR %in% c(2009:2022)) |>
  mutate(EXPCATCHWT = ifelse(is.na(EXPCATCHWT), 0, EXPCATCHWT))

sp_complete_summary <- complete |>
  group_by(DECDEG_BEGLAT, DECDEG_BEGLON, AVGDEPTH) |>
  summarise(
    presence = if_else(any(PRESENCE == 1), 1, 0),      # Presence/absence based on species presence
    species_count = n_distinct(SVSPP),               # Number of distinct species
    sample_size = n()) |>                                  # Sample size (number of rows per location)
  ungroup()

# Select location data and rename columns
data_locations <- sp_complete_summary |>
  dplyr::select(AVGDEPTH, DECDEG_BEGLON, DECDEG_BEGLAT ) |>
  rename(longitude = DECDEG_BEGLON, latitude = DECDEG_BEGLAT, depth = AVGDEPTH)


# Convert data frame to sf object
obs_locations <- st_as_sf(data_locations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# View the sf object
print(obs_locations)
summary(data_locations)

###SIMULATED LOCATIONS###
# Load and process simulated data
survdat_fixed <- readRDS(here("R", "scup", "data", "surv_dat_fixed.rds"))
sim_locations <- survdat_fixed$setdet %>%
  select(x, y, depth, cell) %>%
  mutate(x = x * 1000, y = y * 1000) # Convert from km to m

# Convert to sf object and transform coordinates to WGS84
utm_crs <- 32618  # UTM Zone 18N
simulated_sf <- st_as_sf(sim_locations, coords = c("x", "y"), crs = utm_crs) %>%
  st_transform(crs = 4326)

# Extract coordinates and drop geometry
simulated_locations <- simulated_sf %>%
  mutate(longitude = st_coordinates(.)[,1], latitude = st_coordinates(.)[,2]) #%>%
# st_drop_geometry()

