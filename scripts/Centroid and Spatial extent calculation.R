## FERREIRA, Vitor

## CALCULATION OF SPATIAL EXTENT ##

## POLYGON SYSTEMS ##

## Packages
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

## Importing coordinates
df <- data.frame(
  Longitude = c(
    -76.0067,
    -75.9569,
    -75.9578,
    -75.8763,
    -76.0114,
    -76.0635,
    -75.9613,
    -75.9750,
    -75.8899,
    -76.0197
  ),
  Latitude = c(
    42.0125,
    42.0755,
    42.2013,
    42.0647,
    42.1581,
    42.0815,
    42.0919,
    42.1126,
    42.0992,
    42.1102
  )
)

## Sites sf in lon/lat
pts <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

## EPSG:6933 = World Cylindrical Equal Area
pts_proj <- st_transform(pts, 6933)

## Convex hull, area and centroid
hull <- st_convex_hull(st_union(pts_proj))
area_km2 <- as.numeric(st_area(hull)) / 1e6
centroid_proj <- st_centroid(hull)

## Hull and centroid to lon/lat
hull_ll <- st_transform(hull, 4326)
centroid <- st_transform(centroid_proj, 4326)
centroid_coords <- st_coordinates(centroid)

## Creating map limits
bbox <- st_bbox(hull_ll)
x_range <- bbox["xmax"] - bbox["xmin"]
y_range <- bbox["ymax"] - bbox["ymin"]
x_pad <- max(as.numeric(x_range) * 0.25, 0.5)
y_pad <- max(as.numeric(y_range) * 0.25, 0.5)
xlim_auto <- c(bbox["xmin"] - x_pad, bbox["xmax"] + x_pad)
ylim_auto <- c(bbox["ymin"] - y_pad, bbox["ymax"] + y_pad)

## World map
world <- ne_countries(scale = "medium", returnclass = "sf")

## Plot
map <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = hull_ll, fill = "lightblue", color = "blue",
          alpha = 0.4, linewidth = 1) +
  geom_sf(data = pts, color = "red", size = 2) +
  geom_sf(data = centroid, color = "black", size = 3,
          shape = 4, stroke = 1.5) +
  coord_sf(xlim = xlim_auto, ylim = ylim_auto) +
  labs(subtitle = paste0("Area = ", round(area_km2, 0), " km²")) +
  theme_minimal()

map

## Results
area_km2
centroid_coords
centroid

## LINEAR SYSTEMS ##

## Importing coordinates
df <- data.frame(
  Longitude = c(92.71746, 93.06716, 93.51, 94.18124, 94.52206, 94.63975),
  Latitude = c(29.09403, 29.04379, 29.16218, 29.19571, 29.44924, 29.46494)
)

## Points sf in lon/lat
pts <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

## Defining UTM projection automatically
mean_lon <- mean(df$Longitude)
mean_lat <- mean(df$Latitude)
utm_zone <- floor((mean_lon + 180) / 6) + 1
epsg_utm <- ifelse(mean_lat >= 0, 32600 + utm_zone, 32700 + utm_zone)

## To UTM
pts_proj <- st_transform(pts, epsg_utm)

## Creating line following sites order
# The sites need to be in spatial order of the river
line_proj <- pts_proj %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

## Cumulative length and midpoint
length_km <- as.numeric(st_length(line_proj)) / 1000
midpoint_proj <- st_line_sample(line_proj, sample = 0.5) %>% st_cast("POINT")
midpoint <- st_transform(midpoint_proj, 4326)
midpoint_coords <- st_coordinates(midpoint)

## Back to lon/lat
line_ll <- st_transform(line_proj, 4326)

## Map limits
bbox <- st_bbox(line_ll)
x_range <- bbox["xmax"] - bbox["xmin"]
y_range <- bbox["ymax"] - bbox["ymin"]
x_pad <- max(as.numeric(x_range) * 0.25, 0.5)
y_pad <- max(as.numeric(y_range) * 0.25, 0.5)
xlim_auto <- c(bbox["xmin"] - x_pad, bbox["xmax"] + x_pad)
ylim_auto <- c(bbox["ymin"] - y_pad, bbox["ymax"] + y_pad)

## Map
world <- ne_countries(scale = "medium", returnclass = "sf")

map <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(data = line_ll, color = "blue", linewidth = 1) +
  geom_sf(data = pts, color = "red", size = 2) +
  geom_sf(data = midpoint, color = "black", size = 3,
          shape = 4, stroke = 1.5) +
  coord_sf(xlim = xlim_auto, ylim = ylim_auto) +
  labs(subtitle = paste0("Cumulative length = ",
                         round(length_km, 1), " km")) +
  theme_minimal()

map

## Results
length_km
midpoint_coords
midpoint

## OVERLAP BETWEEN REFERENCE AND IMPACTED SYSTEMS ##

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

## Choose "linear" or "polygon"
system_type <- "polygon"

## Coordinates
df_ref <- data.frame(
  Longitude = c(
    -53.59863, -53.55191, -53.55061, -53.54802
  ),
  Latitude = c(
    -22.86000, -22.84625, -22.82250, -22.79375
  )
)

df_imp <- data.frame(
  Longitude = c(
    -53.32740, -53.23527, -53.21191, -53.19634
  ),
  Latitude = c(
    -22.76375, -22.74250, -22.74000, -22.70750
  )
)

stopifnot(system_type %in% c("linear", "polygon"))

## Projection: UTM for linear systems; EPSG:6933 for polygons
if (system_type == "linear") {
  mean_lon <- mean(c(df_ref$Longitude, df_imp$Longitude))
  mean_lat <- mean(c(df_ref$Latitude, df_imp$Latitude))
  utm_zone <- floor((mean_lon + 180) / 6) + 1
  epsg_proj <- ifelse(mean_lat >= 0, 32600 + utm_zone, 32700 + utm_zone)
} else {
  epsg_proj <- 6933
}

## Create points and spatial geometry
create_system <- function(df) {
  pts <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)
  pts_proj <- st_transform(pts, epsg_proj)
  
  geometry <- if (system_type == "polygon") {
    st_convex_hull(st_union(pts_proj))
  } else {
    st_geometry(
      pts_proj %>%
        summarise(do_union = FALSE) %>%
        st_cast("LINESTRING")
    )
  }
  
  list(pts = pts, geometry = geometry)
}

ref <- create_system(df_ref)
imp <- create_system(df_imp)

## Geometric intersection
overlap <- st_intersects(
  ref$geometry,
  imp$geometry,
  sparse = FALSE
)[1, 1]

## Back to longitude/latitude
ref_ll <- st_transform(ref$geometry, 4326)
imp_ll <- st_transform(imp$geometry, 4326)

## Automatic map limits
bbox <- st_bbox(st_union(ref_ll, imp_ll))
x_pad <- max(as.numeric(bbox["xmax"] - bbox["xmin"]) * 0.25, 0.5)
y_pad <- max(as.numeric(bbox["ymax"] - bbox["ymin"]) * 0.25, 0.5)

xlim_auto <- c(bbox["xmin"] - x_pad, bbox["xmax"] + x_pad)
ylim_auto <- c(bbox["ymin"] - y_pad, bbox["ymax"] + y_pad)

## Map
world <- ne_countries(scale = "medium", returnclass = "sf")

map <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_sf(
    data = ref_ll, color = "blue",
    fill = if (system_type == "polygon") "lightblue" else NA,
    alpha = if (system_type == "polygon") 0.4 else 1,
    linewidth = 1
  ) +
  geom_sf(
    data = imp_ll, color = "red",
    fill = if (system_type == "polygon") "lightcoral" else NA,
    alpha = if (system_type == "polygon") 0.4 else 1,
    linewidth = 1
  ) +
  geom_sf(data = ref$pts, color = "blue", size = 2) +
  geom_sf(data = imp$pts, color = "red", size = 2) +
  coord_sf(xlim = xlim_auto, ylim = ylim_auto) +
  labs(
    subtitle = paste0(
      tools::toTitleCase(system_type),
      " system | Blue = Reference; Red = Impacted | Overlap = ",
      ifelse(overlap, "YES", "NO")
    )
  ) +
  theme_minimal()

map
overlap

## END ##