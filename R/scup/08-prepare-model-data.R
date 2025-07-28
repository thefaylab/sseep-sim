
library(sf)



# Directories
sseep.analysis <- "C:/Users/croman1/Desktop/UMassD/sseep-analysis"
dist.dat <- here("data", "rds", "dists")
survdat <- here("data", "rds", "survdat")
surv.prods <- here("data", "rds", "surv-prods")
mods.data <- here("data", "rds", "surv-prods", "mods_data", "scup")
plots <- here("outputs", "plots")



# Parameters
species   <- "scup"
season    <- "fall"
ages      <- 0:7
years     <- 1:15
nsims     <- 1:100
nsurveys  <- 25
ids     <- sprintf("%03d", nsims)
chunk_size <- 20
chunks <- split(nsims, ceiling(nsims / chunk_size))


survdat_precl <- map(ids, ~readRDS(here(survdat, sprintf("%s_%s_%s_25_precl_survey.rds", species, season, .x))))

# area weights for each strata
strata_wts <- readRDS(here(sseep.analysis, "data", "rds", "active_strata_wts.rds")) |>
  rename(strat = STRATUM)


add_latlon <- function(df, crs_proj = 32618) {
  df$x_m <- df$x * 1000  # convert from km to meters
  df$y_m <- df$y * 1000
  sf_obj <- st_as_sf(df, coords = c("x_m", "y_m"), crs = crs_proj, remove = FALSE)
  sf_obj <- st_transform(sf_obj, crs = 4326)
  coords <- st_coordinates(sf_obj)
  df$lon <- coords[, 1]
  df$lat <- coords[, 2]

  return(df)
}



# survdat_precl is your original list of 100 population data.tables
survdat_precl <- lapply(survdat_precl, function(df) {
  df <- as.data.frame(df)     # in case it's a data.table
  df <- add_latlon(df)        # append lon/lat
  return(df)
})



tb_precl <- map2_dfr(survdat_precl, seq_along(survdat_precl), function(surv, pop_num) {
  surv |> #it is already a data.table
    as_tibble() |>
    filter(strat %in% unlist(strat)) |>
    group_by(sim, year, strat) |>
    left_join(strata_wts, by = "strat")|>
    mutate(scenario = "Preclusion",
           pop = pop_num,
           SEASON = "FALL") |>
    rename(YEAR = year,
           X = x,
           Y = y,
           DECDEG_LAT = lat,
           DECDEG_LON = lon,
           AVGDEPTH = depth,
           STRATUM = strat) |>
    select(set,sim,YEAR,STRATUM,AVGDEPTH,X,Y,DECDEG_LAT,DECDEG_LON,cell,AREA_CODE,N,scenario,pop,SEASON)
})




#plot data
tb_precl1 <- tb_precl |> filter(pop==1, sim==1)
ggplot(tb_precl1, aes(x = X, y = Y)) +
  geom_point(aes(color = N), size = 1.8) +
  scale_color_viridis_c(option = "plasma", trans = "log1p") +  # log transform for better visual spread
  coord_fixed() +
  labs(title = "Abundance per Tow (N)",
       subtitle = "Preclusion Scenario",
       color = "Abundance (N)") +
  theme_minimal()


##LOOP to filter by pop and sim

for (p in unique(tb_precl$pop)) {
for (s in unique(tb_precl$sim)) {

    # Subset data
    sub <- tb_precl %>%
      filter(pop == p, sim == s)

    file_name <- sprintf("precl_pop%03d_sim%02d.rds", p, s)
    file_path <- file.path(mods.data, file_name) # Save full path
    write_rds(sub, file_path, compress = "gz")  # Write to disk immediately
  }
}



