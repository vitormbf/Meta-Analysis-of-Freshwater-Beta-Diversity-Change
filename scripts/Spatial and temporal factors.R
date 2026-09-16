## FERREIRA, Vitor

## CONFOUNDING FACTORS ##

## PACKAGES
library(readxl)
library(writexl)
library(dplyr)
library(purrr)
library(geosphere)
library(lme4)
library(lmerTest)

data <- read_excel("Extraction dataset.xlsx",
                   na = c("", "NA", "NaN", "N/A", "#N/A"))

original_cols_spatial <- names(data)

metrics_spatial <- c(tax_beta_Mean = "Taxonomic",
                     fun_beta_Mean = "Functional",
                     phy_beta_Mean = "Phylogenetic")

design_col_spatial <- if ("Design" %in% names(data)) "Design" else "Deisign"
num_clean_spatial <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

data_spatial <- data %>%
  mutate(Row_ID = row_number(),
         lnRR_Group = trimws(as.character(lnRR_Group)),
         across(all_of(c(names(metrics_spatial), "Longitude", "Latitude")),
                num_clean_spatial))

## LNRR PAIRS ##

calculate_lnRR_spatial <- function(metric_col, metric_name) {
  treatments <- data_spatial %>%
    filter(lnRR_Group %in% c("Treatment", "Impacted"),
           .data[[metric_col]] > 0) %>%
    mutate(Metric = metric_name,
           Treatment_Row_ID = Row_ID,
           Treatment_Impact_degree = Impact_degree,
           Treatment_Region = Region,
           Treatment_Longitude = Longitude,
           Treatment_Latitude = Latitude,
           Treatment_Geometry = Geometry_Polygonkm2_Linearkm,
           Treatment_Beta = .data[[metric_col]])
  
  controls <- data_spatial %>%
    filter(lnRR_Group %in% c("Control", "Reference"),
           .data[[metric_col]] > 0) %>%
    transmute(`Study case_ID`,
              Metric = metric_name,
              Control_Row_ID = Row_ID,
              Control_Impact_degree = Impact_degree,
              Control_Region = Region,
              Control_Longitude = Longitude,
              Control_Latitude = Latitude,
              Control_Geometry = Geometry_Polygonkm2_Linearkm,
              Control_Beta = .data[[metric_col]])
  
  treatments %>%
    inner_join(controls,
               by = c("Study case_ID", "Metric"),
               relationship = "many-to-many") %>%
    mutate(lnRR = log(Treatment_Beta / Control_Beta))
}

lnRR_pairs_spatial <- map2_dfr(names(metrics_spatial),
                               metrics_spatial,
                               calculate_lnRR_spatial) %>%
  select(all_of(original_cols_spatial), Metric,
         Treatment_Row_ID, Treatment_Impact_degree, Treatment_Region,
         Treatment_Longitude, Treatment_Latitude, Treatment_Geometry,
         Treatment_Beta, Control_Row_ID, Control_Impact_degree,
         Control_Region, Control_Longitude, Control_Latitude,
         Control_Geometry, Control_Beta, lnRR) %>%
  arrange(`Study case_ID`, Metric, Treatment_Row_ID, Control_Row_ID)

write_xlsx(lnRR_pairs_spatial, "lnRR_all_pairs_spatial.xlsx")

## COMMON VARIABLES ##

analysis_data_spatial <- lnRR_pairs_spatial %>%
  mutate(Design_clean = recode(trimws(as.character(.data[[design_col_spatial]])),
                               "Interspersed" = "Spatially interspersed",
                               "Spatiotemporal" = "Spatiotemporal",
                               "Segregated" = "Spatially segregated",
                               .default = trimws(as.character(.data[[design_col_spatial]]))),
         Geometry_clean = trimws(as.character(Geometry_Polygonkm2_Linearkm)),
         lnRR_value = as.numeric(lnRR),
         spatial_extent_value = num_clean_spatial(Spatial_extent_Final),
         temporal_window_value = num_clean_spatial(Before_After_Temporal_window))

## GEOGRAPHIC DISTANCE ##

pair_distance_data_spatial <- analysis_data_spatial %>%
  filter(Design_clean == "Spatially segregated",
         !is.na(Treatment_Longitude), !is.na(Treatment_Latitude),
         !is.na(Control_Longitude), !is.na(Control_Latitude)) %>%
  mutate(pair_distance_km = geosphere::distHaversine(
    cbind(Control_Longitude, Control_Latitude),
    cbind(Treatment_Longitude, Treatment_Latitude),
    r = 6371000
  ) / 1000) %>%
  filter(pair_distance_km > 0)

## SPATIAL EXTENT ##

spatial_extent_data <- analysis_data_spatial %>%
  filter(Geometry_clean == "Polygon",
         Design_clean %in% c("Spatially interspersed", "Spatiotemporal"),
         spatial_extent_value > 0)

## ABSOLUTE LATITUDE ##

latitude_data <- analysis_data_spatial %>%
  mutate(paired_latitude = rowMeans(
    cbind(Treatment_Latitude, Control_Latitude), na.rm = TRUE),
    paired_latitude = replace(
      paired_latitude, is.nan(paired_latitude), NA_real_),
    latitude_value = abs(if_else(
      Design_clean == "Spatially segregated",
      paired_latitude, Latitude))) %>%
  filter(!is.na(latitude_value))
 
## TEMPORAL WINDOW ##

temporal_window_data <- analysis_data_spatial %>%
  filter(Design_clean == "Spatiotemporal",
         temporal_window_value > 0)

## MIXED MODELS ##

m_lat <- lmer(lnRR_value ~ latitude_value + (1 | `Study case_ID`),
              data = latitude_data)

m_extent <- lmer(lnRR_value ~ log(spatial_extent_value) +
                   (1 | `Study case_ID`),
                 data = spatial_extent_data)

m_distance <- lmer(lnRR_value ~ log(pair_distance_km) +
                     (1 | `Study case_ID`),
                   data = pair_distance_data_spatial)

m_time <- lmer(lnRR_value ~ temporal_window_value +
                 (1 | `Study case_ID`),
               data = temporal_window_data)

summary(m_lat)
confint(m_lat)

summary(m_extent)
confint(m_extent)

summary(m_distance)
confint(m_distance)

summary(m_time)
confint(m_time)

## NUMBER OF ARTICLES / STUDY CASES / lnRR IN EACH MODEL

## LATITUDE
latitude_data %>%
  filter(
    !is.na(lnRR_value),
    !is.na(latitude_value),
    !is.na(`Study case_ID`)
  ) %>%
  summarise(
    n_Articles = n_distinct(Article_ID),
    n_Study_cases = n_distinct(`Study case_ID`),
    n_lnRR = n()
  )

## SPATIAL EXTENT
spatial_extent_data %>%
  filter(
    !is.na(lnRR_value),
    !is.na(spatial_extent_value),
    !is.na(`Study case_ID`)
  ) %>%
  summarise(
    n_Articles = n_distinct(Article_ID),
    n_Study_cases = n_distinct(`Study case_ID`),
    n_lnRR = n()
  )

## GEOGRAPHIC DISTANCE
pair_distance_data_spatial %>%
  filter(
    !is.na(lnRR_value),
    !is.na(pair_distance_km),
    !is.na(`Study case_ID`)
  ) %>%
  summarise(
    n_Articles = n_distinct(Article_ID),
    n_Study_cases = n_distinct(`Study case_ID`),
    n_lnRR = n()
  )

## TEMPORAL WINDOW
temporal_window_data %>%
  filter(
    !is.na(lnRR_value),
    !is.na(temporal_window_value),
    !is.na(`Study case_ID`)
  ) %>%
  summarise(
    n_Articles = n_distinct(Article_ID),
    n_Study_cases = n_distinct(`Study case_ID`),
    n_lnRR = n()
  )

## NUMBER OF OBSERVATIONS AND UNIQUE MODERATOR VALUES

library(dplyr)

count_model_values <- function(model, variable_name) {
  
  mf <- model.frame(model)
  
  data.frame(
    Variable = variable_name,
    n_observations = nrow(mf),
    n_unique_values = length(unique(mf[[2]]))
  )
}

unique_values_table <- bind_rows(
  
  count_model_values(
    m_lat,
    "Absolute latitude"
  ),
  
  count_model_values(
    m_extent,
    "Spatial extent"
  ),
  
  count_model_values(
    m_distance,
    "Geographic distance"
  ),
  
  count_model_values(
    m_time,
    "Temporal window"
  )
)

unique_values_table

## END ##
