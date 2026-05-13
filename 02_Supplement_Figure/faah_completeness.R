# Load packages
library(tidyverse)
library(ggplot2)
library(scales) 

# Import data
data <- data.frame(
  Threshold = c("com70", "com75", "com80", "com85", "com90", "com95",
                "total_70", "total_75", "total_80", "total_85", "total_90", "total_95"),
  F.longum = c(2130, 2023, 1870, 1640, 1247, 606,
               2614, 2480, 2272, 1951, 1451, 680),
  F.prausnitzii = c(4, 3, 3, 2, 1, 1,
                    1776, 1678, 1559, 1410, 1150, 556),
  F.duncaniae = c(7, 7, 5, 3, 2, 1,
                  1448, 1343, 1208, 1085, 904, 479),
  F.faecis = c(4, 4, 3, 2, 2, 1,
               1332, 1226, 1119, 1004, 858, 523),
  F.un1 = c(2, 2, 1, 0, 0, 0,
            502, 420, 360, 271, 190, 59),
  F.un2 = c(0, 0, 0, 0, 0, 0,
            391, 363, 339, 311, 267, 140),
  F.hattorii = c(0, 0, 0, 0, 0, 0,
                 222, 217, 208, 198, 166, 83),
  F.un3 = c(195, 150, 110, 84, 65, 34,
            208, 159, 116, 90, 68, 35),
  F.butyricigenerans = c(128, 120, 116, 109, 87, 46,
                         134, 124, 120, 113, 91, 49),
  F.un4 = c(65, 62, 60, 54, 48, 30,
            67, 63, 61, 55, 49, 30),
  F.wellingii = c(0, 0, 0, 0, 0, 0,
                  67, 61, 58, 53, 45, 27),
  F.langellae = c(0, 0, 0, 0, 0, 0,
                  64, 62, 58, 55, 44, 16)
)

completeness_data <- data[1:6, ]
total_data <- data[7:12, ]

# Convert data
completeness_long <- completeness_data %>%
  pivot_longer(
    cols = -Threshold,
    names_to = "Species",
    values_to = "Yes"
  ) %>%
  mutate(
    Species = factor(Species, levels = c(
      "F.longum", "F.prausnitzii", "F.duncaniae", 
      "F.faecis", "F.un1", "F.un2", "F.hattorii", 
      "F.un3", "F.butyricigenerans", "F.un4", 
      "F.wellingii", "F.langellae"
    )),
    Threshold = factor(Threshold, levels = c("com70", "com75", "com80", 
                                             "com85", "com90", "com95")),
    Threshold_num = str_extract(Threshold, "\\d+")
  )

total_long <- total_data %>%
  pivot_longer(
    cols = -Threshold,
    names_to = "Species",
    values_to = "Total"
  ) %>%
  mutate(
    Species = factor(Species, levels = levels(completeness_long$Species)),
    Threshold_num = str_extract(Threshold, "\\d+"),
    Threshold = str_replace(Threshold, "total_", "com")
  )

# Merge data
merged_data <- completeness_long %>%
  left_join(total_long %>% select(Species, Threshold_num, Total), 
            by = c("Species", "Threshold_num")) %>%
  select(-Threshold_num) %>%
  mutate(
    No = Total - Yes,
    Yes_pct = Yes / Total,
    No_pct = No / Total
  )

plot_data <- merged_data %>%
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
    x_pos = (Species_num - 1) * 6 + as.numeric(factor(Threshold, levels = c("com70", "com75", "com80", "com85", "com90", "com95")))
  )

# Color mapping
color_mapping <- c(
  # com70
  "com70_No" = "#EFEDF5",
  "com70_Yes" = "#BCBDDC",
  
  # com75
  "com75_No" = "#EFEDF5",
  "com75_Yes" = "#9E9AC8",
  
  # com80
  "com80_No" = "#EFEDF5",
  "com80_Yes" = "#807DBA",
  
  # com85
  "com85_No" = "#EFEDF5",
  "com85_Yes" = "#6A51A3",
  
  # com90
  "com90_No" = "#EFEDF5",
  "com90_Yes" = "#54278F",
  
  # com95
  "com95_No" = "#EFEDF5",
  "com95_Yes" = "#3F007D"
)

legend_color_mapping <- c(
  "com70" = "#BCBDDC",
  "com75" = "#9E9AC8",
  "com80" = "#807DBA",
  "com85" = "#6A51A3",
  "com90" = "#54278F",
  "com95" = "#3F007D"
)

# Set labels
species_total_labels <- merged_data %>%
  arrange(Species, Threshold) %>%
  group_by(Species) %>%
  summarize(
    total_70 = Total[Threshold == "com70"],
    total_75 = Total[Threshold == "com75"],
    total_80 = Total[Threshold == "com80"],
    total_85 = Total[Threshold == "com85"],
    total_90 = Total[Threshold == "com90"],
    total_95 = Total[Threshold == "com95"],
    .groups = "drop"
  ) %>%
  mutate(
    total_label = paste("n=", total_70, ", ", total_75, ", ", total_80, ", ",
                        total_85, ", ", total_90, ", ", total_95, sep = ""),
    Species_num = as.numeric(Species),
    x_center = (Species_num - 1) * 6 + 3.5
  )

species_count <- length(levels(plot_data$Species))
x_breaks <- seq(from = 3.5, by = 6, length.out = species_count)

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

# Print and save plot
p <- ggplot(plot_data, aes(x = x_pos, y = Percentage, fill = Threshold_Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  labs(
    title = "", 
    x = "Species cluster", 
    y = "Percentage of genomes (COG annotation)", 
    fill = "    FAAH
    Completeness"
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
    data = species_total_labels,
    aes(x = x_center, y = 1.02, label = total_label),
    inherit.aes = FALSE,
    size = 0.75,
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
    labels = c("70%", "75%", "80%", "85%", "90%", "95%"),
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
    xintercept = seq(6.5, by = 6, length.out = species_count - 1),
    linetype = "dashed",
    color = "gray50",
    alpha = 0.7
  )


print(p)
ggsave(filename = "FAAH_completeness_all.pdf", plot = p, width = 8, height = 6)
