#
## FERREIRA ET AL. SCRIPT 4

library(tidyr)

pa_data <- lnRR_pairs %>%
  filter(Community_data == "Presence-absence")

pa_boot <- bootstrap_lnRR(
  data = pa_data,
  moderator = NULL
)

pa_boot$summary

###

quant_existing <- lnRR_pairs %>%
  filter(Community_data != "Presence-absence")

sc_recalc <- data %>%
  filter(
    !is.na(AB_tax_beta_Mean) | !is.na(AB_fun_beta_Mean)
  ) %>%
  distinct(`Study case_ID`) %>%
  pull(`Study case_ID`)

sc_recalc
length(sc_recalc)

quant_recalc <- data %>%
  filter(`Study case_ID` %in% c("S8", "S40", "S54", "S67")) %>%
  mutate(Row_ID = row_number()) %>%
  pivot_longer(
    cols = c(AB_tax_beta_Mean, AB_fun_beta_Mean),
    names_to = "Metric_AB",
    values_to = "Beta_AB"
  ) %>%
  mutate(
    Metric = recode(
      Metric_AB,
      "AB_tax_beta_Mean" = "Taxonomic",
      "AB_fun_beta_Mean" = "Functional"
    )
  ) %>%
  filter(!is.na(Beta_AB), Beta_AB > 0) %>%
  select(-Metric_AB)

quant_recalc <- quant_recalc %>%
  filter(lnRR_Group == "Treatment") %>%
  rename(Treatment_Beta = Beta_AB) %>%
  inner_join(
    quant_recalc %>%
      filter(lnRR_Group == "Control") %>%
      select(`Study case_ID`, Metric, Control_Beta = Beta_AB),
    by = c("Study case_ID", "Metric"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    lnRR = log(Treatment_Beta / Control_Beta),
    Community_data = "Quantitative"
  )

lnRR_pairs_quantitative <- bind_rows(
  quant_existing,
  quant_recalc
)

###

quant_boot <- bootstrap_lnRR(
  data = lnRR_pairs_quantitative,
  moderator = NULL
)

quant_boot$summary


###


community_boot_summary <- bind_rows(
  pa_boot$summary %>%
    mutate(Level = "Presence-absence"),
  quant_boot$summary %>%
    mutate(Level = "Quantitative")
) %>%
  mutate(
    Moderator = "Data type",
    Direction = case_when(
      CI_low > 0 ~ "Differentiation",
      CI_high < 0 ~ "Homogenization",
      TRUE ~ "No clear direction"
    )
  )

community_plot_data <- community_boot_summary %>%
  mutate(
    Level_label = case_when(
      Level == "Presence-absence" ~ paste0(
        Level, "  (",
        n_Articles, " A; ",
        n_Study_cases, " SC; ",
        n_lnRR_total, " lnRR)"
      ),
      TRUE ~ paste0(
        Level, "  (",
        n_Articles, "; ",
        n_Study_cases, "; ",
        n_lnRR_total, ")"
      )
    ),
    Level_label = factor(Level_label, levels = rev(Level_label))
  )

community_colors <- c(
  "Presence-absence" = "#7A3E5C",
  "Quantitative" = "#7A3E5C"
)

community_boot_plot <- ggplot(
  community_plot_data,
  aes(x = Mean_lnRR, y = Level_label)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high,
      color = Level
    ),
    height = 0.3,
    linewidth = 1.2
  ) +
  geom_point(
    aes(fill = Level),
    shape = 21,
    color = "white",
    stroke = 1,
    size = 4
  ) +
  facet_grid(
    Moderator ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_manual(values = community_colors) +
  scale_fill_manual(values = community_colors) +
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
    panel.grid.minor = element_blank()
  )

community_boot_plot

ggsave(
  filename = "Community data sensitivity.tiff",
  plot = community_boot_plot,
  width = 9,
  height = 2,
  dpi = 600,
  device = ragg::agg_tiff,
  compression = "lzw"
)



