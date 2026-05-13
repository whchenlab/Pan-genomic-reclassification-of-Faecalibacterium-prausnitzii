# Load packages
library(tidyverse)
library(ggplot2)
library(scales)

# Import data
data <- data.frame(
  Threshold = c("90 threshold", "80 threshold", "70 threshold", "60 threshold", "50 threshold", "COG", "total number"),
  F.longum = c(2203, 2213, 2224, 2215, 2252, 2130, 2614),
  F.prausnitzii = c(4, 4, 4, 49, 54, 4, 1776),
  F.duncaniae = c(1, 2, 3, 3, 9, 7, 1448),
  F.faecis = c(1, 1, 2, 2, 5, 4, 1332),
  F.un1 = c(2, 2, 2, 2, 3, 2, 502),
  F.un2 = c(0, 0, 0, 0, 3, 0, 391),
  F.hattorii = c(0, 0, 0, 0, 4, 0, 222),
  F.un3 = c(0, 0, 194, 194, 195, 195, 208),
  F.butyricigenerans = c(0, 128, 128, 128, 128, 128, 134),
  F.un4 = c(0, 0, 0, 61, 63, 65, 67),
  F.wellingii = c(0, 0, 0, 0, 0, 0, 67),
  F.langellae = c(0, 0, 0, 0, 0, 0, 64)
)
total_counts <- data[data$Threshold == "total number", -1]
data <- data[data$Threshold != "total number", ]

# Convert data
data_long <- data %>%
  pivot_longer(
    cols = -Threshold,
    names_to = "Species",
    values_to = "Yes"
  ) %>%
  mutate(
    # set species order
    Species = factor(Species, levels = c(
      "F.longum", "F.prausnitzii", "F.duncaniae", 
      "F.faecis", "F.un1", "F.un2", "F.hattorii", 
      "F.un3", "F.butyricigenerans", "F.un4", 
      "F.wellingii", "F.langellae"
    )),
    Threshold = factor(Threshold, levels = c("90 threshold", "80 threshold", "70 threshold", 
                                             "60 threshold", "50 threshold", "COG"))
  )

# Add total number data
data_long <- data_long %>%
  rowwise() %>%
  mutate(
    Total = total_counts[[Species]],
    No = Total - Yes,
    # 计算百分比
    Yes_pct = Yes / Total,
    No_pct = No / Total
  ) %>%
  ungroup()

# Create plot data
plot_data <- data_long %>%
  select(Species, Threshold, Yes, No, Yes_pct, No_pct, Total) %>%
  pivot_longer(
    cols = c(Yes, No),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  pivot_longer(
    cols = c(Yes_pct, No_pct),
    names_to = "Category_pct",
    values_to = "Percentage"
  ) %>%
  filter(
    (Category == "Yes" & Category_pct == "Yes_pct") |
      (Category == "No" & Category_pct == "No_pct")
  ) %>%
  mutate(
    Category = factor(Category, levels = c("No", "Yes")),
    Percentage_label = ifelse(Count > 0 & Percentage > 0, 
                              percent(Percentage, accuracy = 0.1), 
                              "")
  )

# Set species group
plot_data <- plot_data %>%
  mutate(
    Threshold_Category = paste(Threshold, Category, sep = "_"),
    Species_num = as.numeric(Species),
    Threshold_num = as.numeric(Threshold),
    x_pos = (Species_num - 1) * 7 + Threshold_num
  )

# Color mapping
color_mapping <- c(
  # 90 threshold
  "90 threshold_No" = "#E5F5E0",
  "90 threshold_Yes" = "#00441B",
  
  # 80 threshold
  "80 threshold_No" = "#E5F5E0",
  "80 threshold_Yes" = "#006D2C",
  
  # 70 threshold
  "70 threshold_No" = "#E5F5E0",
  "70 threshold_Yes" = "#238B45",
  
  # 60 threshold
  "60 threshold_No" = "#E5F5E0",
  "60 threshold_Yes" = "#41AB5D",
  
  # 50 threshold
  "50 threshold_No" = "#E5F5E0",
  "50 threshold_Yes" = "#74C476",
  
  # COG
  "COG_No" = "#E5F5E0",
  "COG_Yes" = "#99CC33"
)

legend_color_mapping <- c(
  "90 threshold" = "#00441B",
  "80 threshold" = "#006D2C",
  "70 threshold" = "#238B45",
  "60 threshold" = "#41AB5D",
  "50 threshold" = "#74C476",
  "COG" = "#99CC33"
)

# Set labels
species_totals <- data_long %>%
  distinct(Species, Total) %>%
  arrange(Species)
species_count <- length(levels(data_long$Species))
x_breaks <- seq(3.5, by = 7, length.out = species_count)

total_labels_data <- data.frame(
  x = x_breaks,
  label = paste0("n=", species_totals$Total),
  y = 1.02
)

species_labels <- c(
  expression(italic("F. longum")),
  expression(italic("F. prausnitzii")),
  expression(italic("F. duncaniae")),
  expression(italic("F. faecis")),
  expression(italic("F. un1")),
  expression(italic("F. un2")),
  expression(italic("F. hattorii")),
  expression(italic("F. un3")),
  expression(italic("F. butyricigenerans")),
  expression(italic("F. un4")),
  expression(italic("F. wellingii")),
  expression(italic("F. langellae"))
)

# Create plot
ggplot(plot_data, aes(x = x_pos, y = Percentage, fill = Threshold_Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  labs(
    title = "", 
    x = "Species cluster (70 CheckM completeness)", 
    y = "Percentage of genomes", 
    fill = "    FAAH
    Identity"
  ) +
  geom_text(
    aes(label = Percentage_label),
    position = position_stack(vjust = 0.5),
    size = 2,
    color = "white",
    fontface = "bold",
    angle = 90,
    na.rm = TRUE
  ) +
  geom_text(
    data = total_labels_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    vjust = 0
  ) +
  scale_fill_manual(
    values = color_mapping,
    guide = guide_legend(
      override.aes = list(
        fill = legend_color_mapping
      )
    ),
    labels = c("90%", "80%", "70%", "60%", "50%", "COG"),
    breaks = paste(names(legend_color_mapping), "Yes", sep = "_")
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = species_labels,
    limits = c(0.5, max(plot_data$x_pos) + 0.5),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "", size = 12, face = 'plain'),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      family = '',
      size = 9,
      face = 'italic'
    ),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.position = 'right',
    plot.margin = margin(20, 10, 10, 10)
  ) +
  geom_vline(
    xintercept = seq(6.5, by = 7, length.out = species_count - 1),
    linetype = "dashed",
    color = "gray50",
    alpha = 0.7
  )
ggsave(filename = "FAAH_blasti_all.pdf", width = 8, height = 6)
