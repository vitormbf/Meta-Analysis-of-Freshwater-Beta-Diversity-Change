## FERREIRA, Vitor ##

## CLIMATIC HETEROGENEITY ##

library(terra)
library(sf)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(lme4)
library(lmerTest)
library(grid)

terraOptions(progress = 1)

wc_folder <- "C:/Users/Vitor/Desktop/Global Change Biology/Data/WorldClim"
target_res <- "30s"

## PREVIOUSLY EXCLUDED ##

excluded_env <- c("S40", "S97")

bio_vars <- c(
  "bio1",   # Annual Mean Temperature
  "bio4",   # Temperature Seasonality
  "bio5",   # Max Temperature of Warmest Month
  "bio6",   # Min Temperature of Coldest Month
  "bio12",  # Annual Precipitation
  "bio15",  # Precipitation Seasonality
  "bio16",  # Precipitation of Wettest Quarter
  "bio17"   # Precipitation of Driest Quarter
)

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

## SPATIAL EXTENT DATA ##

required_cols <- c(
  "Study case_ID",
  "Design_clean",
  "Longitude",
  "Latitude",
  "spatial_extent_value",
  "lnRR_value"
)

missing_cols <- setdiff(
  required_cols,
  names(spatial_extent_data)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing columns in spatial_extent_data: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

spatial_extent_clean <- spatial_extent_data %>%
  filter(
    Design_clean %in% c(
      "Spatially interspersed",
      "Spatiotemporal"
    ),
    !`Study case_ID` %in% excluded_env
  ) %>%
  mutate(
    Longitude = as.numeric(
      gsub(",", ".", as.character(Longitude))
    ),
    Latitude = as.numeric(
      gsub(",", ".", as.character(Latitude))
    ),
    spatial_extent_value = as.numeric(
      gsub(",", ".", as.character(spatial_extent_value))
    ),
    lnRR_value = as.numeric(
      gsub(",", ".", as.character(lnRR_value))
    )
  )

## CONFIRM THE RETAINED DESIGN AND SAMPLE SIZE ##

table(
  spatial_extent_clean$Design_clean,
  useNA = "ifany"
)

cat(
  "Study cases retained:",
  n_distinct(spatial_extent_clean$`Study case_ID`),
  "\n"
)

cat(
  "lnRR values retained:",
  nrow(spatial_extent_clean),
  "\n"
)

## LOADING WORLDCLIM VARIABLES ##

all_tifs <- list.files(
  wc_folder,
  pattern = "\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

length(all_tifs)
head(all_tifs)

bio_files <- all_tifs[
  str_detect(
    tolower(basename(all_tifs)),
    "bio[_-]?\\d+\\.tif$"
  ) &
    str_detect(
      tolower(all_tifs),
      target_res
    )
]

## FALLBACK: IF THE FOLDER CONTAINS ONLY 30S FILES BUT FILENAMES/PATH DO NOT INCLUDE "30S" ##

if (length(bio_files) == 0) {
  bio_files <- all_tifs[
    str_detect(
      tolower(basename(all_tifs)),
      "bio[_-]?\\d+\\.tif$"
    )
  ]
}

if (length(bio_files) == 0) {
  stop(
    paste0(
      "No WorldClim bioclimatic .tif files found. ",
      "Check wc_folder or unzip the files."
    )
  )
}

bio_num <- as.numeric(
  str_match(
    tolower(basename(bio_files)),
    "bio[_-]?(\\d+)\\.tif$"
  )[, 2]
)

bio_files <- bio_files[order(bio_num)]
bio_num <- bio_num[order(bio_num)]

bio <- rast(bio_files)

names(bio) <- paste0("bio", bio_num)

missing_bio <- setdiff(
  bio_vars,
  names(bio)
)

if (length(missing_bio) > 0) {
  stop(
    paste(
      "Missing bioclimatic variables:",
      paste(missing_bio, collapse = ", ")
    )
  )
}

bio_sel <- bio[[bio_vars]]

bio_sel

plot(bio_sel[[1]])

## EQUIVALENT-AREA RADIUS ##

env_centers <- spatial_extent_clean %>%
  group_by(`Study case_ID`) %>%
  summarise(
    Longitude = first_non_na(Longitude),
    Latitude = first_non_na(Latitude),
    spatial_extent_value = first_non_na(
      spatial_extent_value
    ),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(Longitude),
    !is.na(Latitude),
    !is.na(spatial_extent_value),
    spatial_extent_value > 0
  ) %>%
  mutate(
    radius_km = sqrt(
      spatial_extent_value / pi
    ),
    radius_m = radius_km * 1000
  )

nrow(env_centers)

setdiff(
  spatial_extent_clean %>%
    distinct(`Study case_ID`) %>%
    pull(`Study case_ID`),
  env_centers$`Study case_ID`
)

## CREATE EQUIVALENT AREA BUFFERS ##

make_buffer_aeqd <- function(
    id,
    lon,
    lat,
    area_km2,
    radius_m
) {
  pt <- st_sfc(
    st_point(c(lon, lat)),
    crs = 4326
  )
  
  local_crs <- paste0(
    "+proj=aeqd ",
    "+lat_0=", lat, " ",
    "+lon_0=", lon, " ",
    "+datum=WGS84 +units=m +no_defs"
  )
  
  buffer_ll <- pt %>%
    st_transform(local_crs) %>%
    st_buffer(dist = radius_m) %>%
    st_transform(4326)
  
  st_sf(
    `Study case_ID` = id,
    Longitude = lon,
    Latitude = lat,
    spatial_extent_value = area_km2,
    radius_km = radius_m / 1000,
    geometry = buffer_ll
  )
}

env_buffers <- pmap_dfr(
  list(
    id = env_centers$`Study case_ID`,
    lon = env_centers$Longitude,
    lat = env_centers$Latitude,
    area_km2 = env_centers$spatial_extent_value,
    radius_m = env_centers$radius_m
  ),
  make_buffer_aeqd
)

nrow(env_buffers)

## CHECK BUFFER AREAS ##

buffer_area_check <- env_buffers %>%
  st_transform(6933) %>%
  mutate(
    buffer_area_km2 =
      as.numeric(st_area(geometry)) / 1e6,
    area_difference_km2 =
      buffer_area_km2 - spatial_extent_value,
    area_difference_percent =
      100 *
      area_difference_km2 /
      spatial_extent_value
  ) %>%
  st_drop_geometry()

buffer_area_check

## EXTRACT WORLDCLIM PIXELS INSIDE BUFFERS ##

buffers_vect <- terra::vect(env_buffers)

bio_pixels <- terra::extract(
  bio_sel,
  buffers_vect,
  cells = TRUE,
  touches = TRUE
) %>%
  as_tibble()

id_lookup <- env_buffers %>%
  st_drop_geometry() %>%
  mutate(
    ID = row_number()
  ) %>%
  select(
    ID,
    `Study case_ID`,
    Longitude,
    Latitude,
    spatial_extent_value,
    radius_km
  )

bio_pixels <- bio_pixels %>%
  left_join(
    id_lookup,
    by = "ID"
  ) %>%
  relocate(
    `Study case_ID`,
    ID,
    cell,
    Longitude,
    Latitude,
    spatial_extent_value,
    radius_km
  )

pixel_count <- bio_pixels %>%
  group_by(`Study case_ID`) %>%
  summarise(
    n_pixels = n(),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_pixels)
  )

pixel_count

setdiff(
  env_centers$`Study case_ID`,
  pixel_count$`Study case_ID`
)

## CREATE CLEAN PIXEL MATRIX ##

bio_pixels_matrix <- bio_pixels %>%
  select(
    `Study case_ID`,
    cell,
    all_of(bio_vars)
  ) %>%
  filter(
    if_all(
      all_of(bio_vars),
      ~ !is.na(.x)
    )
  )

## PCA AND ENVIRONMENTAL HETEROGENEITY ##

bio_pixels_clean <- bio_pixels_matrix %>%
  filter(
    if_all(
      all_of(bio_vars),
      ~ !is.na(.x)
    )
  )

pca_clim <- prcomp(
  bio_pixels_clean %>%
    select(all_of(bio_vars)),
  center = TRUE,
  scale. = TRUE
)

var_exp <- pca_clim$sdev^2 /
  sum(pca_clim$sdev^2)

pca_var <- tibble(
  PC = paste0(
    "PC",
    seq_along(pca_clim$sdev)
  ),
  variance_explained = var_exp,
  cumulative_variance = cumsum(var_exp)
)

pca_var

n_pcs <- which(
  pca_var$cumulative_variance >= 0.90
)[1]

pc_cols <- paste0(
  "PC",
  1:n_pcs
)

pca_scores <- predict(
  pca_clim,
  newdata = bio_pixels_clean %>%
    select(all_of(bio_vars))
) %>%
  as_tibble() %>%
  select(
    all_of(pc_cols)
  )

bio_scores <- bind_cols(
  bio_pixels_clean %>%
    select(
      `Study case_ID`,
      cell
    ),
  pca_scores
)

centroids <- bio_scores %>%
  group_by(`Study case_ID`) %>%
  summarise(
    across(
      all_of(pc_cols),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  rename_with(
    ~ paste0(.x, "_centroid"),
    all_of(pc_cols)
  )

bio_scores_centroid <- bio_scores %>%
  left_join(
    centroids,
    by = "Study case_ID"
  )

score_matrix <- as.matrix(
  bio_scores_centroid[, pc_cols]
)

centroid_matrix <- as.matrix(
  bio_scores_centroid[
    ,
    paste0(pc_cols, "_centroid")
  ]
)

bio_scores_centroid$dist_to_centroid <- sqrt(
  rowSums(
    (score_matrix - centroid_matrix)^2
  )
)

climatic_heterogeneity <- bio_scores_centroid %>%
  group_by(`Study case_ID`) %>%
  summarise(
    climatic_heterogeneity =
      mean(
        dist_to_centroid,
        na.rm = TRUE
      ),
    climatic_heterogeneity_sd =
      sd(
        dist_to_centroid,
        na.rm = TRUE
      ),
    climatic_heterogeneity_median =
      median(
        dist_to_centroid,
        na.rm = TRUE
      ),
    n_pixels = n(),
    .groups = "drop"
  ) %>%
  arrange(
    desc(climatic_heterogeneity)
  )

climatic_heterogeneity

## JOIN ENVIRONMENTAL HETEROGENEITY WITH LNRR DATA ##

env_data <- spatial_extent_clean %>%
  left_join(
    climatic_heterogeneity,
    by = "Study case_ID"
  )

env_data %>%
  summarise(
    n_study_cases =
      n_distinct(`Study case_ID`),
    n_study_cases_with_env =
      n_distinct(
        `Study case_ID`[
          !is.na(climatic_heterogeneity)
        ]
      ),
    n_lnRR = n()
  )

## CORRELATION: SPATIAL EXTENT VS ENVIRONMENTAL HETEROGENEITY ##

env_sc_level <- env_data %>%
  distinct(
    `Study case_ID`,
    spatial_extent_value,
    climatic_heterogeneity,
    n_pixels
  ) %>%
  filter(
    !is.na(spatial_extent_value),
    !is.na(climatic_heterogeneity),
    spatial_extent_value > 0
  )

cor.test(
  log10(
    env_sc_level$spatial_extent_value
  ),
  env_sc_level$climatic_heterogeneity
)

plot_extent_het <- ggplot(
  env_sc_level,
  aes(
    x = log10(spatial_extent_value),
    y = climatic_heterogeneity
  )
) +
  geom_point(
    size = 2.3,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  labs(
    x = expression(
      log[10]~"spatial extent (km²)"
    ),
    y = "Climatic heterogeneity"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    )
  )

plot_extent_het

ggsave(
  paste0(
    "worldclim_30s_spatial_extent_vs_",
    "climatic_heterogeneity_",
    "spatially_interspersed.tiff"
  ),
  plot_extent_het,
  width = 5.5,
  height = 4.2,
  units = "in"
)

## MIXED MODEL: ENVIRONMENTAL HETEROGENEITY AND LNRR ##

env_model_data <- env_data %>%
  mutate(
    Study_case =
      as.factor(`Study case_ID`),
    lnRR_value =
      as.numeric(lnRR_value),
    climatic_heterogeneity_z =
      as.numeric(
        scale(climatic_heterogeneity)
      )
  ) %>%
  filter(
    !is.na(Study_case),
    !is.na(lnRR_value),
    !is.na(climatic_heterogeneity_z),
    is.finite(lnRR_value),
    is.finite(climatic_heterogeneity_z)
  )

model_clim_het <- lmer(
  lnRR_value ~
    climatic_heterogeneity_z +
    (1 | Study_case),
  data = env_model_data,
  REML = TRUE
)

summary(model_clim_het)

confint(
  model_clim_het,
  method = "Wald"
)

## OPTIONAL: SAME SUBSET, SPATIAL EXTENT MODEL ##

env_model_data <- env_model_data %>%
  mutate(
    log_spatial_extent_z =
      as.numeric(
        scale(
          log10(spatial_extent_value)
        )
      )
  )

model_extent_same_subset <- lmer(
  lnRR_value ~
    log_spatial_extent_z +
    (1 | Study_case),
  data = env_model_data,
  REML = TRUE
)

summary(model_extent_same_subset)

confint(
  model_extent_same_subset,
  method = "Wald"
)

## PCA BIPLOT ##

set.seed(123)

scores_plot_base <- bio_scores %>%
  select(
    `Study case_ID`,
    PC1,
    PC2
  )

n_sample <- min(
  nrow(scores_plot_base),
  100000
)

scores_plot <- scores_plot_base %>%
  slice_sample(
    n = n_sample
  )

study_centroids_plot <- bio_scores %>%
  group_by(`Study case_ID`) %>%
  summarise(
    PC1 = mean(
      PC1,
      na.rm = TRUE
    ),
    PC2 = mean(
      PC2,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    climatic_heterogeneity,
    by = "Study case_ID"
  )

loadings_plot <- as_tibble(
  pca_clim$rotation[, 1:2],
  rownames = "Variable"
) %>%
  mutate(
    Variable_label = recode(
      Variable,
      "bio1" =
        "Annual mean temp.",
      "bio4" =
        "Temp. seasonality",
      "bio5" =
        "Max temp. warmest month",
      "bio6" =
        "Min temp. coldest month",
      "bio12" =
        "Annual precip.",
      "bio15" =
        "Precip. seasonality",
      "bio16" =
        "Precip. wettest quarter",
      "bio17" =
        "Precip. driest quarter",
      .default = Variable
    )
  )

arrow_scale <- min(
  diff(
    range(
      scores_plot$PC1,
      na.rm = TRUE
    )
  ),
  diff(
    range(
      scores_plot$PC2,
      na.rm = TRUE
    )
  )
) * 0.30

loadings_plot <- loadings_plot %>%
  mutate(
    PC1_arrow =
      PC1 * arrow_scale,
    PC2_arrow =
      PC2 * arrow_scale
  )

pc1_lab <- paste0(
  "PC1 (",
  round(
    pca_var$variance_explained[1] * 100,
    1
  ),
  "%)"
)

pc2_lab <- paste0(
  "PC2 (",
  round(
    pca_var$variance_explained[2] * 100,
    1
  ),
  "%)"
)

pca_biplot <- ggplot() +
  geom_point(
    data = scores_plot,
    aes(
      x = PC1,
      y = PC2
    ),
    color = "grey75",
    alpha = 0.08,
    size = 0.25
  ) +
  geom_point(
    data = study_centroids_plot,
    aes(
      x = PC1,
      y = PC2,
      color = climatic_heterogeneity
    ),
    size = 3,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = loadings_plot,
    aes(
      x = PC1_arrow,
      y = PC2_arrow,
      label = Variable_label
    ),
    color = "firebrick",
    size = 3.2,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  scale_color_viridis_c(
    option = "C",
    name = "Environmental heterogeneity"
  ) +
  labs(
    x = pc1_lab,
    y = pc2_lab
  ) +
  coord_equal() +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 9
    )
  )

pca_biplot

ggsave(
  paste0(
    "worldclim_30s_pca_biplot_",
    "climatic_heterogeneity_",
    "spatially_interspersed.tiff"
  ),
  pca_biplot,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

## PCA PLOT WITH CONVEX HULLS BY STUDY CASE ##

pca_hulls <- scores_plot_base %>%
  group_by(`Study case_ID`) %>%
  filter(n() >= 3) %>%
  slice(
    chull(
      PC1,
      PC2
    )
  ) %>%
  ungroup() %>%
  left_join(
    climatic_heterogeneity,
    by = "Study case_ID"
  )

pca_biplot_hulls <- ggplot() +
  geom_polygon(
    data = pca_hulls,
    aes(
      x = PC1,
      y = PC2,
      group = `Study case_ID`,
      fill = climatic_heterogeneity
    ),
    alpha = 0.25,
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_point(
    data = study_centroids_plot,
    aes(
      x = PC1,
      y = PC2,
      color = climatic_heterogeneity
    ),
    size = 1
  ) +
  scale_fill_viridis_c(
    option = "C",
    guide = "none"
  ) +
  scale_color_viridis_c(
    option = "C",
    name = "Climatic heterogeneity"
  ) +
  labs(
    x = pc1_lab,
    y = pc2_lab
  ) +
  coord_equal() +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 9
    )
  )

pca_biplot_hulls

ggsave(
  "PCA_spatially_interspersed.tiff",
  pca_biplot_hulls,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

## ENVIRONMENTAL HETEROGENEITY ##

env_model_data %>%
  summarise(
    n_Articles = n_distinct(Article_ID, na.rm = TRUE),
    n_Study_cases = n_distinct(`Study case_ID`, na.rm = TRUE),
    n_lnRR = n()
  )