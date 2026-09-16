## FERREIRA, Vitor

## BOOTSTRAP RESAMPLING PROCEDURE ##

## PACKAGES ##
library(readxl)
library(dplyr)
library(purrr)
library(writexl)
library(ragg)
library(ggplot2)

## CALCULATION LNRR FOR EVERY POSSIBLE PAIR WITHIN STUDY CASES ##

data <- read_excel("Extraction dataset.xlsx", na = c("", "NA", "NaN", "N/A", "#N/A"))

original_cols <- names(data)
metrics <- c(tax_beta_Mean = "Taxonomic", fun_beta_Mean = "Functional", phy_beta_Mean = "Phylogenetic")

data <- data %>%
  mutate(Row_ID = row_number(),
         lnRR_Group = trimws(as.character(lnRR_Group)),
         across(all_of(names(metrics)), ~ suppressWarnings(as.numeric(gsub(",", ".", .x)))))

calculate_lnRR <- function(metric_col, metric_name) {
  treatments <- data %>%
    filter(lnRR_Group == "Treatment", !is.na(.data[[metric_col]])) %>%
    mutate(Metric = metric_name,
           Treatment_Row_ID = Row_ID,
           Treatment_Impact_degree = Impact_degree,
           Treatment_Beta = .data[[metric_col]])
  
  controls <- data %>%
    filter(lnRR_Group == "Control", !is.na(.data[[metric_col]])) %>%
    transmute(`Study case_ID`,
              Control_Row_ID = Row_ID,
              Control_Impact_degree = Impact_degree,
              Control_Beta = .data[[metric_col]])
  
  treatments %>%
    inner_join(controls, by = "Study case_ID", relationship = "many-to-many") %>%
    mutate(lnRR = log(Treatment_Beta / Control_Beta))
}

lnRR_pairs <- map2_dfr(names(metrics), metrics, calculate_lnRR) %>%
  select(all_of(original_cols), Metric,
         Treatment_Row_ID, Treatment_Impact_degree, Treatment_Beta,
         Control_Row_ID, Control_Impact_degree, Control_Beta, lnRR) %>%
  arrange(`Study case_ID`, Metric, Treatment_Row_ID, Control_Row_ID)

glimpse(lnRR_pairs)

## OVERALL BOOTSTRAP RESAMPLING PROCEDURE ##

bootstrap_lnRR <- function(data, moderator = NULL, n_iter = 10000, seed = 123,
                           id_col = "Study case_ID", article_col = "Article_ID",
                           effect_col = "lnRR") {
  set.seed(seed)
  
  moderator_name <- if (is.null(moderator)) "Overall" else moderator
  
  bootstrap_data <- data %>%
    transmute(Study_case = .data[[id_col]],
              Article_ID = .data[[article_col]],
              Effect_size = as.numeric(.data[[effect_col]]),
              Moderator = if (is.null(moderator)) "Overall" else as.character(.data[[moderator]])) %>%
    filter(!is.na(Study_case), !is.na(Effect_size), !is.na(Moderator))
  
  results <- lapply(unique(bootstrap_data$Moderator), function(level_name) {
    level_data <- bootstrap_data %>%
      filter(Moderator == level_name)
    
    bootstrap_means <- replicate(n_iter, {
      sampled_data <- level_data %>%
        group_by(Study_case) %>%
        slice_sample(n = 1) %>%
        ungroup()
      
      mean(sampled_data$Effect_size)
    })
    
    summary <- tibble(Moderator = moderator_name,
                      Level = level_name,
                      Mean_lnRR = mean(bootstrap_means),
                      Median_lnRR = median(bootstrap_means),
                      CI_low = as.numeric(quantile(bootstrap_means, 0.025)),
                      CI_high = as.numeric(quantile(bootstrap_means, 0.975)),
                      n_iter = n_iter,
                      n_Articles = n_distinct(level_data$Article_ID, na.rm = TRUE),
                      n_Study_cases = n_distinct(level_data$Study_case),
                      n_lnRR_total = nrow(level_data))
    
    iterations <- tibble(Moderator = moderator_name,
                         Level = level_name,
                         Iteration = seq_len(n_iter),
                         Mean_lnRR = bootstrap_means)
    
    list(summary = summary, iterations = iterations)
  })
  
  list(summary = bind_rows(lapply(results, `[[`, "summary")),
       iterations = bind_rows(lapply(results, `[[`, "iterations")))
}

## RUNNING BOOTSTRAPS ##

overall_boot <- bootstrap_lnRR(data = lnRR_pairs, moderator = NULL, n_iter = 10000)
overall_boot$summary

metric_boot <- bootstrap_lnRR(data = lnRR_pairs, moderator = "Metric", n_iter = 10000)
metric_boot$summary

impact_boot <- bootstrap_lnRR(data = lnRR_pairs, moderator = "Impact", n_iter = 10000)
impact_boot$summary

design_boot <- bootstrap_lnRR(data = lnRR_pairs, moderator = "Design", n_iter = 10000)
design_boot$summary

group_boot <- bootstrap_lnRR(data = lnRR_pairs, moderator = "Group", n_iter = 10000)
group_boot$summary

## PLOTTING MAIN RESULTS ##

all_boot_summary <- bind_rows(
  overall_boot$summary,
  metric_boot$summary,
  impact_boot$summary,
  design_boot$summary,
  group_boot$summary
) %>%
  filter(
    n_Articles >= 10,
    n_Study_cases >= 10,
    n_lnRR_total >= 10
  ) %>%
  mutate(
    Moderator = recode(
      Moderator,
      "Overall" = "Global",
      "Metric" = "Facet",
      "Impact" = "Impact",
      "Design" = "Design",
      "Group" = "Group"
    ),
    Level = recode(
      Level,
      "Overall" = "Global",
      "Segregated" = "Spatially segregated",
      "Interspersed" = "Spatially interspersed", 
      "Spatiotemporal" = "Spatiotemporal"
    ),
    Direction = case_when(
      CI_low > 0 ~ "Differentiation",
      CI_high < 0 ~ "Homogenization",
      TRUE ~ "No clear direction"
    )
  )

section_order <- c(
  "Global",
  "Facet",
  "Impact",
  "Design",
  "Group"
)

level_order <- c(
  "Global",
  "Taxonomic", "Functional",
  "Hydrological alteration", "Pollution/nutrients", "Land use", "Multiple", "Invasion",
  "Spatially segregated", "Spatially interspersed", "Spatiotemporal",
  "Fish", "Invertebrates", "Zooplankton"
)

moderator_colors <- c(
  "Global" = "#9E9E9E",
  "Facet"  = "#3B6EA8",
  "Impact" = "#B65A3A",
  "Design" = "#4E8F74",
  "Group"  = "#8E6AAE"
)

plot_data <- all_boot_summary %>%
  mutate(
    Moderator = factor(Moderator, levels = section_order),
    order_level = match(Level, level_order),
    order_level = if_else(is.na(order_level), 999L, order_level),
    Level_label = if_else(
      Level == "Global",
      paste0(
        Level, "  (",
        n_Articles, " A; ",
        n_Study_cases, " SC; ",
        n_lnRR_total, " lnRR)"
      ),
      paste0(
        Level, "  (",
        n_Articles, "; ",
        n_Study_cases, "; ",
        n_lnRR_total, ")"
      )
    ),
    Plot_ID = paste(Moderator, Level_label, sep = "___")
  ) %>%
  arrange(Moderator, order_level, Level) %>%
  mutate(
    Plot_ID = factor(Plot_ID, levels = rev(unique(Plot_ID)))
  )

all_boot_plot <- ggplot(plot_data, aes(x = Mean_lnRR, y = Plot_ID)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      color = Moderator
    ),
    height = 0.3,
    linewidth = 1.2
  ) +
  geom_point(
    aes(fill = Moderator),
    shape = 21,
    color = "white",
    stroke = 1,
    size = 4
  ) +
  facet_grid(
    Moderator ~ .,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      Moderator = function(x) ifelse(x == "Global", "", x)
    )
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_color_manual(values = moderator_colors) +
  scale_fill_manual(values = moderator_colors) +
  guides(
    color = "none",
    fill = "none"
  ) +
  labs(
    x = "Effect size (lnRR)",
    y = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.right = element_text(
      angle = 90,
      hjust = 0.5, 
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    panel.spacing.y = unit(0.8, "lines")
  )

all_boot_plot

ggsave(
  filename = "Main results.tiff",
  plot = all_boot_plot,
  width = 9,
  height = 6,
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

## PLOTTING LOW-SAMPLE RESULTS ##

low_sample_boot_summary <- bind_rows(
  overall_boot$summary,
  metric_boot$summary,
  impact_boot$summary,
  design_boot$summary,
  group_boot$summary
) %>%
  filter(
    n_Articles <= 9 | n_Study_cases <= 9 | n_lnRR_total <= 9
  ) %>%
  mutate(
    Moderator = recode(
      Moderator,
      "Overall" = "Global",
      "Metric" = "Facet",
      "Impact" = "Impact",
      "Design" = "Design",
      "Group" = "Group"
    ),
    Level = recode(
      Level,
      "Overall" = "Global",
      "Segregated" = "Spatially segregated",
      "Interspersed" = "Spatially interspersed", 
      "Spatiotemporal" = "Spatiotemporal"
    ),
    Direction = case_when(
      CI_low > 0 ~ "Differentiation",
      CI_high < 0 ~ "Homogenization",
      TRUE ~ "No clear direction"
    )
  )

low_sample_level_order <- c(
  "Global",
  "Taxonomic", "Functional", "Phylogenetic",
  "Hydrological alteration", "Pollution/nutrients", "Land use", "Multiple", "Invasion",
  "Spatially segregated", "Spatially interspersed", "Spatiotemporal",
  "Fish", "Invertebrates", "Zooplankton", "Algae", "Phytoplankton", "Macrophyte"
)

low_sample_plot_data <- low_sample_boot_summary %>%
  mutate(
    Moderator = factor(Moderator, levels = section_order),
    order_level = match(Level, low_sample_level_order),
    order_level = if_else(is.na(order_level), 999L, order_level)
  ) %>%
  arrange(Moderator, order_level, Level) %>%
  mutate(
    Level_label = if_else(
      row_number() == 1,
      paste0(
        Level, "  (",
        n_Articles, " A; ",
        n_Study_cases, " SC; ",
        n_lnRR_total, " lnRR)"
      ),
      paste0(
        Level, "  (",
        n_Articles, "; ",
        n_Study_cases, "; ",
        n_lnRR_total, ")"
      )
    ),
    Plot_ID = paste(Moderator, Level_label, sep = "___"),
    Plot_ID = factor(Plot_ID, levels = rev(unique(Plot_ID)))
  )

low_sample_boot_plot <- ggplot(low_sample_plot_data, aes(x = Mean_lnRR, y = Plot_ID)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      color = Moderator
    ),
    height = 0.3,
    linewidth = 1.2
  ) +
  geom_point(
    aes(fill = Moderator),
    shape = 21,
    color = "white",
    stroke = 1,
    size = 4
  ) +
  facet_grid(
    Moderator ~ .,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      Moderator = function(x) ifelse(x == "Global", "", x)
    )
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_color_manual(values = moderator_colors) +
  scale_fill_manual(values = moderator_colors) +
  guides(
    color = "none",
    fill = "none"
  ) +
  labs(
    x = "Effect size (lnRR)",
    y = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.right = element_text(
      angle = 90,
      hjust = 0.5, 
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    panel.spacing.y = unit(0.8, "lines")
  )

low_sample_boot_plot

ggsave(
  filename = "Low sample results.tiff",
  plot = low_sample_boot_plot,
  width = 9,
  height = 7,
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

## GLOBAL BOOTSTRAP WITHOUT FISH AND INVERTEBRATES ##

other_groups_data <- lnRR_pairs %>%
  filter(
    !Group %in% c("Fish", "Invertebrates")
  )

global_without_fish_inv_boot <- bootstrap_lnRR(
  data = other_groups_data,
  moderator = NULL,
  n_iter = 10000
)

global_without_fish_inv_boot$summary

## LEAVE-ONE-OUT BOOTSTRAP ##

leave_one_out_bootstrap_lnRR <- function(data,
                                         n_iter = 10000,
                                         seed = 123,
                                         id_col = "Study case_ID",
                                         article_col = "Article_ID",
                                         effect_col = "lnRR") {
  
  if (!effect_col %in% names(data) && "lnrr" %in% names(data)) {
    data <- data %>% rename(lnRR = lnrr)
    effect_col <- "lnRR"
  }
  
  if (!id_col %in% names(data)) {
    stop(paste0("Column '", id_col, "' not found in data."))
  }
  
  if (!article_col %in% names(data)) {
    stop(paste0("Column '", article_col, "' not found in data."))
  }
  
  if (!effect_col %in% names(data)) {
    stop(paste0("Column '", effect_col, "' not found in data."))
  }
  
  # Same basic cleaning used by the bootstrap
  valid_data <- data %>%
    mutate(
      Study_case = as.character(.data[[id_col]]),
      Article_ID_clean = .data[[article_col]],
      Effect_size = as.numeric(.data[[effect_col]])
    ) %>%
    filter(
      !is.na(Study_case),
      !is.na(Effect_size)
    )
  
  study_cases <- sort(unique(valid_data$Study_case))
  
  # Full-data bootstrap
  full_boot <- bootstrap_lnRR(
    data = data,
    moderator = NULL,
    n_iter = n_iter,
    seed = seed,
    id_col = id_col,
    article_col = article_col,
    effect_col = effect_col
  )
  
  full_summary <- full_boot$summary %>%
    mutate(Analysis = "Full data")
  
  full_mean <- full_summary$Mean_lnRR[1]
  full_ci_low <- full_summary$CI_low[1]
  full_ci_high <- full_summary$CI_high[1]
  
  # Information on what is removed each time
  removed_info <- valid_data %>%
    group_by(Study_case) %>%
    summarise(
      Removed_n_lnRR = n(),
      Removed_n_Articles = n_distinct(Article_ID_clean),
      .groups = "drop"
    )
  
  # Leave-one-study-case-out bootstrap
  loo_summary <- purrr::map_dfr(seq_along(study_cases), function(i) {
    
    sc <- study_cases[i]
    
    data_loo <- data %>%
      filter(as.character(.data[[id_col]]) != sc)
    
    boot_loo <- bootstrap_lnRR(
      data = data_loo,
      moderator = NULL,
      n_iter = n_iter,
      seed = seed + i,
      id_col = id_col,
      article_col = article_col,
      effect_col = effect_col
    )
    
    boot_loo$summary %>%
      mutate(
        Removed_Study_case = sc,
        Analysis = "Leave-one-out"
      )
  }) %>%
    left_join(
      removed_info,
      by = c("Removed_Study_case" = "Study_case")
    ) %>%
    mutate(
      Full_Mean_lnRR = full_mean,
      Full_CI_low = full_ci_low,
      Full_CI_high = full_ci_high,
      Delta_Mean_lnRR = Mean_lnRR - full_mean,
      Abs_Delta_Mean_lnRR = abs(Delta_Mean_lnRR),
      CI_overlaps_zero = CI_low <= 0 & CI_high >= 0,
      Inference = case_when(
        CI_low > 0 ~ "Differentiation",
        CI_high < 0 ~ "Homogenization",
        TRUE ~ "Overlaps zero"
      )
    ) %>%
    arrange(desc(Abs_Delta_Mean_lnRR))
  
  list(
    full = full_summary,
    leave_one_out = loo_summary
  )
}

loo_global <- leave_one_out_bootstrap_lnRR(
  data = lnRR_pairs,
  n_iter = 10000,
  seed = 123,
  id_col = "Study case_ID",
  article_col = "Article_ID",
  effect_col = "lnRR"
)

loo_global$full

loo_global$leave_one_out %>%
  select(
    Removed_Study_case,
    Mean_lnRR,
    CI_low,
    CI_high,
    Delta_Mean_lnRR,
    Abs_Delta_Mean_lnRR,
    Removed_n_lnRR,
    Removed_n_Articles,
    Inference
  ) %>%
  arrange(desc(Abs_Delta_Mean_lnRR)) %>%
  print(n = Inf)

loo_plot_data <- loo_global$leave_one_out %>%
  mutate(
    Removed_Study_case = reorder(Removed_Study_case, Mean_lnRR)
  )

loo <- ggplot(loo_plot_data, aes(x = Removed_Study_case, y = Mean_lnRR)) +
  geom_hline(
    yintercept = loo_global$full$Mean_lnRR[1],
    linetype = "dashed", color = "red"
  ) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high)) +
  coord_flip() +
  labs(
    x = "Removed study case",
    y = "Mean lnRR after removing one study case (95% CI)",
  ) +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )

ggsave(
  filename = "Leave-one-out.tiff",
  plot = loo,
  width = 8,
  height = 12,
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

loo_global$leave_one_out %>%
  count(n_Study_cases)

loo_delta <- loo_global$leave_one_out %>%
  mutate(
    Removed_Study_case = reorder(Removed_Study_case, Delta_Mean_lnRR)
  ) %>%
  ggplot(aes(x = Removed_Study_case, y = Delta_Mean_lnRR)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  coord_flip() +
  labs(
    x = "Removed study case",
    y = "Change in mean lnRR after removal",
  ) +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )

loo_delta

ggsave(
  filename = "Leave-one-out DELTA.tiff",
  plot = loo_delta,
  width = 8,
  height = 12,
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)

## END ##





