# load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Import data
data <- data.frame(
  row.names = c("total_95", "total_90", "total_85", "total_80", "total_75", "total_70"),
  F.longum = c(680, 1451, 1951, 2272, 2480, 2614),
  F.prausnitzii = c(556, 1150, 1410, 1559, 1678, 1776),
  F.duncaniae = c(479, 904, 1085, 1208, 1343, 1448),
  F.faecis = c(523, 858, 1004, 1119, 1226, 1332),
  F.un1 = c(59, 190, 271, 360, 420, 502),
  F.un2 = c(140, 267, 311, 339, 363, 391),
  F.hattorii = c(83, 166, 198, 208, 217, 222),
  F.un3 = c(35, 68, 90, 116, 159, 208),
  F.butyricigenerans = c(49, 91, 113, 120, 124, 134),
  F.un4 = c(30, 49, 55, 61, 63, 67),
  F.wellingii = c(27, 45, 53, 58, 61, 67),
  F.langellae = c(16, 44, 55, 58, 62, 64)
)

# Completeness order
data$Threshold <- c(95, 90, 85, 80, 75, 70)

# Convert data
data_long <- data %>%
  pivot_longer(cols = -Threshold, names_to = "Species", values_to = "Genome_Number")

# Color mapping
color_mapping <- c(
  "F.longum" = '#FF7256',
  "F.prausnitzii" = '#FFA500',
  "F.duncaniae" = "#FDD80A",
  "F.faecis" = "#556B2F",
  "F.hattorii" = "#A1D99B",
  "F.butyricigenerans" = "#7FFFd4",
  "F.wellingii" = "#8EE5EE",
  "F.langellae" = "#6495ED",
  "F.un1" = "#B23AEE",
  "F.un2" = "#8B008B",
  "F.un3" = "#FF1493",
  "F.un4" = "#8B0A50"
)
data_long$Species <- factor(data_long$Species, levels = names(color_mapping))

# Print and save
ggplot(data_long, aes(x = Threshold, y = Genome_Number, color = Species, group = Species)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = color_mapping, name = "Species") +
  scale_x_reverse(breaks = seq(70, 95, by = 5), labels = seq(70, 95, by = 5)) +
  labs(
    x = "CheckM Completeness Threshold",
    y = "Genome Number",
    title = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
  ) +
  guides(color = guide_legend(ncol = 1))

ggsave("cluster_number.pdf", width = 6, height = 5)
