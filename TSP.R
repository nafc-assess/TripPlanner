
library(TSP)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Rstrap)

## Random selection of a trip --------------------------------------------------

all_sets <- setdet[setdet$rec == 5, ]
all_sets$trip_id <- with(all_sets, paste0(survey.year, "-", vessel, "-", trip))
all_sets <- all_sets[!duplicated(paste0(all_sets$trip_id, "-", all_sets$set)), ]
unique_trips <- unique(all_sets$trip_id)

time <- all_sets$time.mid
all_sets$hour <- as.numeric(substr(time, 1, nchar(time) - 2))
all_sets$minute <- as.numeric(substr(time, nchar(time) - 1, nchar(time)))
all_sets$date <- with(all_sets, ISOdate(year, month, day, hour, minute))
all_sets$long.start <- -all_sets$long.start

for (t in tail(unique_trips, 50)) {

  one_trip <- t # "2023-76-33" # sample(unique_trips, 1)
  sets <- all_sets[all_sets$trip_id == one_trip, ]
  sets <- sets[order(sets$date), ] # make sure to sort by date to re-construct path

  ## Extract coordinates and solve TSP -------------------------------------------

  coords <- data.frame(lon = sets$long.start, lat = sets$lat.start)
  coords_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  coords_utm <- st_transform(coords_sf, crs = 32622)

  actual_path <- data.frame(st_coordinates(coords_utm))
  dist_mat <- as.matrix(dist(actual_path))
  tsp <- TSP(dist_mat)
  actual_solution <- solve_TSP(tsp, method = "identity")
  tsp_solution <- solve_TSP(tsp)
  tsp_order <- as.integer(tsp_solution)
  tsp_path <- actual_path[tsp_order, ]

  actual_tot_dist <- attr(actual_solution, "tour_length")
  tsp_tot_dist <- attr(tsp_solution, "tour_length")

  ## Map paths -------------------------------------------------------------------

  world_map <- ne_coastline(scale = "large", returnclass = "sf")
  bbox <- c(xmin = -75, xmax = -30, ymin = 30, ymax = 65) # NL bbox
  nl_map <- st_crop(world_map, bbox) |>
    st_transform(crs = 32622)

  plot(coords_utm, type = "o", pch = 16, main = one_trip)
  plot(st_geometry(nl_map), add = TRUE, col = "lightgrey")
  plot(st_geometry(ne_countries(country = "canada")), add = TRUE)
  with(tsp_path, lines(X, Y, col = "red"))
  legend("topright", lty = c(1, 1), col = c("black", "red"),
         legend = c(paste0("Actual (", round(actual_tot_dist / 1000), " km)"),
                    paste0("TSP (", round(tsp_tot_dist / 1000), " km)")))
  box()

}

## Check improvement across trips ----------------------------------------------

comp_dist <- function(trip) {
  one_trip <- sample(unique_trips, 1)
  sets <- all_sets[all_sets$trip_id == trip, ]
  sets <- sets[order(sets$date), ] # make sure to sort by date to re-construct path

  coords <- data.frame(lon = sets$long.start, lat = sets$lat.start)
  coords_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
  coords_utm <- st_transform(coords_sf, crs = 32622)

  actual_path <- data.frame(st_coordinates(coords_utm))
  dist_mat <- as.matrix(dist(actual_path))
  tsp <- TSP(dist_mat)
  actual_solution <- solve_TSP(tsp, method = "identity")
  tsp_solution <- solve_TSP(tsp)
  tsp_order <- as.integer(tsp_solution)
  tsp_path <- actual_path[tsp_order, ]

  actual_tot_dist <- attr(actual_solution, "tour_length")
  tsp_tot_dist <- attr(tsp_solution, "tour_length")
  data.frame(actual = actual_tot_dist, tsp = tsp_tot_dist)
}

dists <- lapply(unique_trips, comp_dist)
dists <- do.call(rbind, dists)

hist((dists$actual - dists$tsp) / 1000, breaks = 50, xlab = "Difference in distance traveled (km)",
     main = "Potential reduction in travel distance using a TSP algorithim")
box()
