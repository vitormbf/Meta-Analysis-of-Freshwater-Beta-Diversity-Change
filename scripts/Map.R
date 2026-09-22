library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

clean_num <- function(x) as.numeric(gsub("[^0-9.-]", "", as.character(x)))

## CLEAN COORDINATES ##

map_data <- data %>%
  mutate(
    Longitude = clean_num(Longitude),
    Latitude = clean_num(Latitude)
  ) %>%
  filter(
    !is.na(Longitude),
    !is.na(Latitude),
    abs(Longitude) <= 180,
    abs(Latitude) <= 90
  ) %>%
  distinct(Longitude, Latitude)

## CONVERT POINTS TO SF ##

pts <- st_as_sf(
  map_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)


## WORLD MAP ##

world <- ne_countries(scale = "medium", returnclass = "sf")

graticule <- st_graticule(
  lon = seq(-180, 180, 60),
  lat = seq(-60, 60, 30),
  crs = 4326
)

lon_labels <- data.frame(
  lon = seq(-180, 180, 60),
  lat = -88,
  label = c("180°", "120°W", "60°W", "0°", "60°E", "120°E", "180°")
)

lat_labels <- data.frame(
  lon = -180,
  lat = seq(-60, 60, 30),
  label = c("60°S", "30°S", "0°", "30°N", "60°N")
)

## MAP ##

map <- ggplot() +
  geom_sf(
    data = graticule,
    color = "grey80",
    linewidth = 0.25
  ) +
  geom_sf(
    data = world,
    fill = "grey75",
    color = "grey75",
    linewidth = 0.2
  ) +
  geom_sf(
    data = pts,
    shape = 21,
    fill = "#B23A48",
    color = "white",
    size = 2.4,
    stroke = 0.5,
    alpha = 0.85
  ) +
  geom_text(
    data = lon_labels,
    aes(x = lon, y = lat, label = label),
    size = 3.2,
    color = "grey35",
    vjust = 2.4
  ) +
  geom_text(
    data = lat_labels,
    aes(x = lon, y = lat, label = label),
    size = 3.2,
    color = "grey35",
    hjust = 1.35
  ) +
  coord_sf(
    crs = "+proj=eqearth",
    default_crs = st_crs(4326),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 15, r = 20, b = 35, l = 45)
  )

map

ggsave(
  "map.tiff",
  plot = map,
  width = 7,
  height = 5,
  dpi = 600,
  compression = "lzw"
)



map_data %>%
  summarise(n_plotted_points = n())


data %>%
  summarise(
    n_sampled_regions = n_distinct(Region)
  )

## END ##
