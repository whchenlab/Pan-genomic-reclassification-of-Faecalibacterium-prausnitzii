library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# Function: process data
process_data <- function(file_path, species_name, desired_genes, title) {

  data <- read.csv(file_path, header = TRUE)
  

  gene_labels <- ifelse(grepl("^[a-zA-Z]", desired_genes), desired_genes, "")
  
  # Extract strain name
  strain_names <- as.character(data$菌株)
  if(any(duplicated(strain_names))) {
    strain_names <- make.unique(strain_names, sep = "_")
  }
  
  # Extract gene data and status
  gene_cols <- 2:(length(desired_genes) + 1)
  gene_data <- data[, gene_cols]
  colnames(gene_data) <- desired_genes
  status_col <- as.character(data$Status)
  
  # Create plot data
  plot_data <- gene_data %>%
    mutate(
      Strain = strain_names,
      Status = status_col,
      RowID = 1:n(),
      Species = species_name
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
  
  return(plot_data)
}

# Function: create plot
create_heatmap <- function(plot_data, title, show_x_axis = TRUE, show_y_axis = TRUE, show_y_title = TRUE) {
  # Strain cluster
  cluster_data <- plot_data %>%
    select(Strain, Gene, IsNull) %>%
    pivot_wider(names_from = Gene, values_from = IsNull) %>%
    as.data.frame()
  
  rownames(cluster_data) <- cluster_data$Strain
  cluster_data <- cluster_data[, -1, drop = FALSE]
  cluster_data <- cluster_data %>% mutate(across(everything(), ~ as.numeric(!.)))
  
  if(nrow(cluster_data) > 1) {
    dist_matrix <- dist(cluster_data, method = "binary")
    hc <- hclust(dist_matrix, method = "ward.D2")
    ordered_strains <- rownames(cluster_data)[hc$order]
  } else {
    ordered_strains <- rownames(cluster_data)
  }
  plot_data$Strain <- factor(plot_data$Strain, levels = ordered_strains)
  
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
  
  gene_labels <- ifelse(grepl("^[a-zA-Z]", levels(plot_data$Gene)), 
                        levels(plot_data$Gene), "")
  
  # create plot
  p <- ggplot(plot_data, aes(x = Gene, y = Strain, fill = Color_Factor)) +
    geom_tile(color = "gray30", linewidth = 0.1) +
    scale_fill_manual(
      values = c(
        "Not Available" = "white",
        "GalNac Cluster Absent" = '#FD2598',
        "Insertion Direction Unknown" = '#1D73FF',
        "GalNac Forward Insertion" = '#7FC97F',
        "GalNac Reverse Insertion" = '#FDD80A'
      ),
      name = "GalNac Cluster Status",
      drop = FALSE
    ) +
    scale_x_discrete(
      breaks = levels(plot_data$Gene),
      labels = gene_labels
    ) +
    labs(
      title = title,
      x = if(show_x_axis) "Target Genes" else "",
      y = if(show_y_title) "Strains" else ""
    ) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6,
        margin = margin(t = 2, b = 2)
      ),
      axis.text.y = element_text(size = 0.5),
      legend.position = "none",  
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  if (!show_x_axis) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  if (!show_y_axis) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  return(p)
}

# Define gene list
genes_v1 <- c(
  "dinB", "2530", "lacC", "2532", "2533", "ptsH", "gatZ-kbaZ", "nagA", 
  "agaS", "agaF", "2539", "gatY-kbaY", "GH109", "rhaR", "agaV", "agaC", 
  "agaD", "pgmB", "2547", "malX", "2549", "2550", "immR_6", "2552", 
  "2553", "2554", "ndoA_2", "2556", "Int-Tn_4", "2558", "2559", "2560"
)

genes_v2 <- c(
  "1666", "RspR", "1668", "1669", "1670", "ptsH", "gatZ-kbaZ", "nagA", 
  "agaS", "agaF", "1676", "gatY-kbaY", "GH109", "rhaR", "agaV", "agaC", 
  "agaD", "1683", "1684", "XerC", "1686", "1687", "1688", "1689", "YwiE", 
  "1691", "1692", "1693", "1694"
)

# Import data
fl_data <- process_data("fl_dsv.csv", "F.hattorii", genes_v2, 
                        "F.hattorii GalNac Cluster Status")
fb_data <- process_data("fb_dsv.csv", "F.wellingii", genes_v2, 
                        "F.wellingii GalNac Cluster Status")
num_fl_strains <- length(unique(fl_data$Strain))
num_fb_strains <- length(unique(fb_data$Strain))
height_ratios <- c(num_fl_strains, num_fb_strains)

# Create heatmap
p_fl <- create_heatmap(fl_data, "F.longum", 
                       show_x_axis = FALSE, show_y_axis = FALSE, show_y_title = TRUE)
p_fb <- create_heatmap(fb_data, "F.butyricigenerans", 
                       show_x_axis = TRUE, show_y_axis = FALSE, show_y_title = TRUE)
combined_heatmaps <- p_fl / p_fb +
  plot_layout(ncol = 1, heights = height_ratios)

# Save as pdf
ggsave("flfb_heatmaps.pdf", combined_heatmaps, 
       width = 4, height = 8, device = "pdf")

