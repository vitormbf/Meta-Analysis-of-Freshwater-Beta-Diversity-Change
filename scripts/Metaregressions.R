## FERREIRA, Vitor

## METAREGRESSIONS ##

## PACKAGES ##
library(readxl)
library(dplyr)
library(tibble)
library(writexl)
library(metafor)
library(broom)

data <- read_excel("Extraction dataset.xlsx", na = c("", "NA", "NaN", "N/A", "#N/A"))

original_cols <- names(data)

metric_info <- tibble(
  Metric = c("Taxonomic", "Functional", "Phylogenetic"),
  mean_col = c("tax_beta_Mean", "fun_beta_Mean", "phy_beta_Mean"),
  sd_col = c("tax_beta_SD", "fun_beta_SD", "phy_beta_SD"),
  n_col = c("tax_beta_N", "fun_beta_N", "phy_beta_N")
)

cols_to_numeric <- unique(c(metric_info$mean_col, metric_info$sd_col, metric_info$n_col))

data <- data %>%
  mutate(
    Row_ID = row_number(),
    lnRR_Group = trimws(as.character(lnRR_Group)),
    across(
      all_of(cols_to_numeric),
      ~ suppressWarnings(as.numeric(gsub(",", ".", as.character(.x))))
    )
  )

make_meta_pairs <- function(metric_name, mean_col, sd_col, n_col) {
  
  treatments <- data %>%
    filter(
      lnRR_Group == "Treatment",
      !is.na(.data[[mean_col]]),
      !is.na(.data[[sd_col]]),
      !is.na(.data[[n_col]]),
      .data[[mean_col]] > 0,
      .data[[sd_col]] >= 0,
      .data[[n_col]] > 0
    ) %>%
    mutate(
      Metric = metric_name,
      Treatment_Row_ID = Row_ID,
      Treatment_Impact_degree = Impact_degree,
      Treatment_Beta = .data[[mean_col]],
      Treatment_SD = .data[[sd_col]],
      Treatment_N = .data[[n_col]]
    )
  
  controls <- data %>%
    filter(
      lnRR_Group == "Control",
      !is.na(.data[[mean_col]]),
      !is.na(.data[[sd_col]]),
      !is.na(.data[[n_col]]),
      .data[[mean_col]] > 0,
      .data[[sd_col]] >= 0,
      .data[[n_col]] > 0
    ) %>%
    transmute(
      `Study case_ID`,
      Metric = metric_name,
      Control_Row_ID = Row_ID,
      Control_Impact_degree = Impact_degree,
      Control_Beta = .data[[mean_col]],
      Control_SD = .data[[sd_col]],
      Control_N = .data[[n_col]]
    )
  
  treatments %>%
    inner_join(
      controls,
      by = c("Study case_ID", "Metric"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      lnRR = log(Treatment_Beta / Control_Beta),
      vi = (Treatment_SD^2 / (Treatment_N * Treatment_Beta^2)) +
        (Control_SD^2 / (Control_N * Control_Beta^2)),
      sei = sqrt(vi)
    )
}

meta_pairs <- lapply(seq_len(nrow(metric_info)), function(i) {
  make_meta_pairs(
    metric_name = metric_info$Metric[i],
    mean_col = metric_info$mean_col[i],
    sd_col = metric_info$sd_col[i],
    n_col = metric_info$n_col[i]
  )
}) %>%
  bind_rows() %>%
  filter(
    !is.na(lnRR),
    !is.na(vi),
    vi > 0
  ) %>%
  mutate(
    Study_case = `Study case_ID`
  )

overall_meta_model <- rma.mv(
  yi = lnRR,
  V = vi,
  random = ~ 1 | Study_case,
  data = meta_pairs,
  method = "REML"
)

summary(overall_meta_model)

## TESTING MODERATORS ##

metric_meta_model <- rma.mv(
  yi = lnRR,
  V = vi,
  mods = ~ 0 + Metric,
  random = ~ 1 | Study_case,
  data = meta_pairs,
  method = "REML"
)

impact_meta_model <- rma.mv(
  yi = lnRR,
  V = vi,
  mods = ~ 0 + Impact,
  random = ~ 1 | Study_case,
  data = meta_pairs,
  method = "REML"
)

design_meta_model <- rma.mv(
  yi = lnRR,
  V = vi,
  mods = ~ 0 + Design,
  random = ~ 1 | Study_case,
  data = meta_pairs,
  method = "REML"
)

group_meta_model <- rma.mv(
  yi = lnRR,
  V = vi,
  mods = ~ 0 + Group,
  random = ~ 1 | Study_case,
  data = meta_pairs,
  method = "REML"
)

summary(metric_meta_model)
summary(impact_meta_model)
summary(design_meta_model)
summary(group_meta_model)

## EGGERS TEST FOR PUBLICATION BIAS ##

library(clubSandwich)

egger_meta_model <- rma.mv(
  yi = lnRR,
  V = vi,
  mods = ~ sei,
  random = ~ 1 | Study_case,
  data = meta_pairs,
  method = "REML"
)

summary(egger_meta_model)

egger_meta_robust <- robust(
  egger_meta_model,
  cluster = meta_pairs$Study_case,
  clubSandwich = TRUE
)

egger_meta_robust

## END ##