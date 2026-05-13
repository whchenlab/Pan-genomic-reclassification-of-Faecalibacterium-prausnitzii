# Prepare packages
library(ggplot2)
library(dplyr)

# 1. Define species name
clusters <- c("F.longum", "F.prausnitzii", "F.duncaniae", "F.faecis", 
              "F.hattorii", "F.butyricigenerans", "F.wellingii", "F.langellae",
              "F.un1", "F.un2", "F.un3", "F.un4")
# 2. Exclude unnamed species
filtered_data <- genome_info_updated %>%
  filter(ANIcluster %in% clusters)

# 3. Color mapping
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
filtered_data$ANIcluster <- factor(filtered_data$ANIcluster, 
                                   levels = names(color_mapping))
missing_clusters <- setdiff(names(color_mapping), unique(filtered_data$ANIcluster))
if(length(missing_clusters) > 0) {
  print(paste("Find absent species cluster:", 
              paste(missing_clusters, collapse = ", ")))
  # Only show present species
  color_mapping <- color_mapping[names(color_mapping) %in% unique(filtered_data$ANIcluster)]
  filtered_data$ANIcluster <- factor(filtered_data$ANIcluster, 
                                     levels = names(color_mapping))
}

# 4. Count species cluster size
cluster_counts <- filtered_data %>%
  group_by(ANIcluster) %>%
  summarise(n = n()) %>%
  ungroup()
filtered_data <- filtered_data %>%
  left_join(cluster_counts, by = "ANIcluster")

# 5. Draw density plot
p <- ggplot(filtered_data, aes(x = Completeness, color = ANIcluster)) +
  geom_density(aes(y = after_stat(count)/n), linewidth = 1) +
  scale_color_manual(values = color_mapping) +
  coord_cartesian(xlim = c(70, 100)) +
  labs(
    x = "Completeness",
    y = "Genome density",
    color = "Species"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

# 6. Save as pdf
print(p)
ggsave("cluster_density_v2.pdf", p, width = 6, height = 5, dpi = 300)
