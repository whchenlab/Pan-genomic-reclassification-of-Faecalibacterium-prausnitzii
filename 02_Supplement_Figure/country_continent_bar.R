# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Process >90 Completeness genomes
genome_info_updated <- genome_info_updated %>%
  filter(Completeness >= 90)

# Sum country count
count_data <- genome_info_updated %>%
  mutate(
    Country_group = ifelse(
      is.na(Country) & is.na(Continent), 
      "Unknown", 
      as.character(Country)
    ),
    Continent_group = ifelse(
      is.na(Country) & is.na(Continent), 
      "N/A", 
      as.character(Continent)
    )
  ) %>%
  group_by(Country_group, Continent_group) %>%
  summarise(Genome_Number = n(), .groups = "drop")

# Set order
count_data <- count_data %>%
  mutate(
    sort_group = ifelse(Country_group == "Unknown", 2, 1),
    sort_value = ifelse(Country_group == "Unknown", 
                        -1,  # Unknown固定值
                        -Genome_Number)  # 其他按数量降序
  ) %>%
  arrange(sort_group, sort_value)

count_data$Country_group <- factor(
  count_data$Country_group,
  levels = rev(unique(count_data$Country_group))
)

# Color mapping
continent_colors <- c(
  "Asia" = "#FF6666",
  "Europe" = "#377EB8",
  "North America" = "#A1D99B",
  "South America" = "#BEAED4",
  "Africa" = "#FDC086",
  "Oceania" = "#FFFF99",
  "N/A" = "#D3D3D3"
)

# Print and save plot
p <- ggplot(count_data, aes(x = Genome_Number, 
                            y = Country_group,
                            fill = Continent_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    name = "Continent",
    values = continent_colors,
    na.value = "gray50"
  ) +
  labs(
    title = "",
    x = "Genome Number",
    y = "Country",
    fill = "Continent"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  geom_text(
    aes(label = Genome_Number), 
    hjust = -0.2, 
    size = 2,
    color = "black"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)))
print(p)
ggsave("country_continent_bar_90completeness.pdf", p, width = 6, height = 4, dpi = 300)
