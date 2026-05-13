# Load required libraries and data
library(ggplot2)
library(dplyr)
matrix <- read.csv("gene_presence_absence.csv")

# Convert to presence/absence matrix
strain_cols <- 4:ncol(matrix)
presence_matrix <- matrix[, strain_cols, drop = FALSE]
presence_matrix[presence_matrix != "" & !is.na(presence_matrix)] <- 1
presence_matrix[presence_matrix == "" | is.na(presence_matrix)] <- 0
presence_matrix <- as.data.frame(lapply(presence_matrix, as.numeric))

# Count number for different type of pan-gene
gene_counts <- rowSums(presence_matrix)
n_genomes <- ncol(presence_matrix)
pangenome_stats <- data.frame(
  category = factor(
    c("core", "softcore", "shell", "cloud"),
    levels = c("core", "softcore", "shell", "cloud")
  ),
  count = c(
    sum(gene_counts >= 0.90 * n_genomes),
    sum(gene_counts >= 0.80 * n_genomes & gene_counts < 0.90 * n_genomes),
    sum(gene_counts >= 0.15 * n_genomes & gene_counts < 0.80 * n_genomes),
    sum(gene_counts < 0.15 * n_genomes)
  ),
  percentage_range = c(
    "90-100%",
    "80-90%",
    "15-80%",
    "0-15%"
  )
) %>%
  mutate(
    percentage = count / sum(count) * 100,
    cum_percentage = cumsum(percentage),
    mid_percentage = cum_percentage - percentage/2,
    legend_text = paste0(
      category, "\n",
      percentage_range, "\n",
      "n=", count
    ),
    label_text = paste0(round(percentage, 1), "%")
  )

# Set position and color for pie plot
pie_x <- 0.5
pie_y <- 0.55
pie_r <- 0.5
pie_colors <- c(
  "#e41a1c",
  "#4daf4a",
  "#377eb8",
  "#984ea3"
)

# Create frequency plot
p <- ggplot(data.frame(gene_counts), aes(x = gene_counts)) +
  geom_histogram(
    binwidth = 1,
    fill = "#f0f0f0",
    color = "black",
    alpha = 0.7
  ) +
  scale_fill_manual(
    values = pie_colors,
    name = "Gene Categories",
    labels = pangenome_stats$legend_text
  ) +
  labs(
    x = "Number of genomes with gene",
    y = "Number of genes",
    title = "Faecalibacterium genus Pangenome"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(1.5, "cm"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8, lineheight = 1.2)
  )

# Calculate actual position for pie plot
x_lim <- layer_scales(p)$x$get_limits()
y_lim <- layer_scales(p)$y$get_limits()
x_length <- x_lim[2] - x_lim[1]
y_length <- y_lim[2] - y_lim[1]
aspect_ratio <- y_length / x_length

pie_x_actual <- x_lim[1] + pie_x * x_length
pie_y_actual <- y_lim[1] + pie_y * y_length
pie_r_actual <- pie_r * x_length

# Create pie plot
pie_plot <- ggplot(pangenome_stats, aes(x = 1, y = percentage, fill = category)) +
  geom_col(width = 1, color = "white", size = 0.7) +
  coord_polar(theta = "y") +
  geom_text(
    aes(
      x = 1.4,
      label = label_text
    ),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  expand_limits(x = c(0.5, 1.5)) +
  scale_fill_manual(values = pie_colors, guide = "none") +
  theme_void()

# Add pie plot and legend
pie_grob <- ggplotGrob(pie_plot)
p <- p + annotation_custom(
  grob = pie_grob,
  xmin = pie_x_actual - pie_r_actual,
  xmax = pie_x_actual + pie_r_actual,
  ymin = pie_y_actual - pie_r_actual * aspect_ratio,
  ymax = pie_y_actual + pie_r_actual * aspect_ratio
)
p <- p + geom_point(
  data = pangenome_stats,
  aes(x = -Inf, y = -Inf, fill = category),
  shape = 21,
  size = 12,
  color = "black",
  show.legend = TRUE
)

# Save as pdf
print(p)
ggsave(
  filename = "pie_freq_fall_v2.pdf",
  width = 6,
  height = 4,
  dpi = 300
)
