# Load required packages and data
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
data <- read.csv("fl_80_with_status.csv", header = TRUE)

# Define gene name, only show galnac-gene
desired_genes <- c(
  "1666", "RspR", "1668", "1669", "1670", "ptsH", "gatZ-kbaZ", "nagA", 
  "agaS", "agaF", "1676", "gatY-kbaY", "GH109", "rhaR", "agaV", "agaC", 
  "agaD", "1683", "1684", "XerC", "1686", "1687", "1688", "1689", "YwiE", 
  "1691", "1692", "1693", "1694"
)
gene_labels <- ifelse(grepl("^[a-zA-Z]", desired_genes), desired_genes, "")

# Extract strain name
strain_names <- as.character(data$菌株)
if(any(duplicated(strain_names))) {
  strain_names <- make.unique(strain_names, sep = "_")
}

# Extract gene presence/absence status and function status
gene_cols <- 2:30
gene_data <- data[, gene_cols]
colnames(gene_data) <- desired_genes
status_col <- as.character(data$Status)

# Create dataframe for heatmap
plot_data <- gene_data %>%
  mutate(
    Strain = strain_names,
    Status = status_col,
    RowID = 1:n()
  ) %>%
  pivot_longer(
    cols = all_of(desired_genes),
    names_to = "Gene",
    values_to = "Value"
  ) %>%
  mutate(
    IsNull = is.na(Value) | Value == "null" | Value == "",
    Color = case_when(
      IsNull ~ "white",
      Status == "Forward" ~ '#7FC97F',
      Status == "Reverse" ~ '#FDD80A',
      Status == "Absent" ~ '#FD2598',
      Status == "Unknown" ~ '#1D73FF',
      TRUE ~ "white"
    )
  )
plot_data$Gene <- factor(plot_data$Gene, levels = desired_genes)

# Strain cluster
cluster_data <- plot_data %>%
  select(Strain, Gene, IsNull) %>%
  pivot_wider(names_from = Gene, values_from = IsNull) %>%
  as.data.frame()

rownames(cluster_data) <- cluster_data$Strain
cluster_data <- cluster_data[, -1, drop = FALSE]
cluster_data <- cluster_data %>% mutate(across(everything(), ~ as.numeric(!.)))

dist_matrix <- dist(cluster_data, method = "binary")
hc <- hclust(dist_matrix, method = "ward.D2")
ordered_strains <- rownames(cluster_data)[hc$order]
plot_data$Strain <- factor(plot_data$Strain, levels = ordered_strains)

# Color mapping for gene and function status
plot_data$Color_Factor <- factor(
  plot_data$Color,
  levels = c("white", '#FD2598', '#1D73FF', '#7FC97F', '#FDD80A'),
  labels = c(
    "Not Available",
    "GalNac Cluster Absent",
    "Insertion Direction Unknown",
    "GalNac Forward Insertion",
    "GalNac Reverse Insertion"
  )
)
color_mapping <- c(
  "Not Available" = "white",
  "GalNac Cluster Absent" = '#FD2598',
  "Insertion Direction Unknown" = '#1D73FF',
  "GalNac Forward Insertion" = '#7FC97F',
  "GalNac Reverse Insertion" = '#FDD80A'
)

# Create heatmap and save as pdf
p <- ggplot(plot_data, aes(x = Gene, y = Strain, fill = Color_Factor)) +
  geom_tile(color = "transparent", linewidth = 0.1) +
  scale_fill_manual(
    values = color_mapping,
    name = "GalNac Cluster Status",
    drop = FALSE
  ) +
  scale_x_discrete(
    breaks = desired_genes,
    labels = gene_labels
  ) +
  labs(
    title = "F.longum GalNac Cluster Status",
    x = "Target Genes",
    y = "Strains"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      margin = margin(t = 2, b = 2)
    ),
    axis.text.y = element_text(size = 0.2),
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(
    title = "GalNac Cluster Status",
    nrow = 2,
    byrow = TRUE
  ))
ggsave("fl_heatmap.pdf", p, 
       width = 6, height =8 , device = "pdf")
print(p)

