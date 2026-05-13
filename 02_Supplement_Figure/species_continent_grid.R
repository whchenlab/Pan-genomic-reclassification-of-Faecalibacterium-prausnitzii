# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Filter >90 Completeness genomes
genome_info_updated <- genome_info_updated %>%
  filter(Completeness >= 90)

# Set order
species_order <- c(
  "F.longum", "F.prausnitzii", "F.duncaniae", "F.faecis", "F.un1", "F.un2",
  "F.hattorii", "F.un3", "F.butyricigenerans", "F.un4", "F.wellingii", "F.langellae"
)
continent_order <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")

# Color mapping
continent_colors <- c(
  "Asia" = "#FF6666",
  "Europe" = "#377EB8",
  "North America" = "#A1D99B",
  "South America" = "#BEAED4",
  "Africa" = "#FDC086",
  "Oceania" = "#FFFF99"
)

# Process data
plot_data <- genome_info_updated %>%
  filter(!is.na(Continent) & !is.na(ANIcluster)) %>%
  filter(ANIcluster %in% species_order) %>%
  group_by(Continent, ANIcluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(Continent = continent_order, ANIcluster = species_order, fill = list(count = 0))

plot_data$ANIcluster <- factor(plot_data$ANIcluster, levels = rev(species_order))
plot_data$Continent <- factor(plot_data$Continent, levels = continent_order)

# Create plot
p <- ggplot(plot_data, aes(x = Continent, y = ANIcluster)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 15)),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  geom_point(aes(size = count, color = Continent), alpha = 0.8) +
  scale_color_manual(values = continent_colors, guide = "none") +
  scale_size_continuous(
    name = "Number of genomes",
    range = c(1, 6),
    breaks = function(x) {
      max_val <- max(x)
      if (max_val <= 5) {
        seq(0, max_val, 1)
      } else if (max_val <= 10) {
        seq(0, max_val, 2)
      } else if (max_val <= 20) {
        seq(0, max_val, 5)
      } else {
        seq(0, max_val, ceiling(max_val/5))
      }
    }
  ) +
  geom_text(
    aes(label = ifelse(count > 0, count, "")),
    size = 1.5,
    color = "black",
    fontface = "bold"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(color = "black"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  labs(title = "")

print(p)
ggsave("species_continent_grid_90completeness.pdf", p, width = 6, height = 4, dpi = 300, bg = "white")
