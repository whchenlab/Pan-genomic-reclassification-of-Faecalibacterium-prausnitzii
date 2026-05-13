# Import
ani_data <- read.table("ani_matrix_genus.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)[,1:3]
colnames(ani_data) <- c("query", "subject", "ANI")
ani_data[, 1:2] <- lapply(ani_data[, 1:2], function(x) sub("\\.fna$", "", x))
# Preprocess
bidirectional_data <- rbind(
  ani_data,
  data.frame(
    query = ani_data$subject,
    subject = ani_data$query,
    ANI = ani_data$ANI
  )
)
# Process inconsistent bidirectional records
library(dplyr)
bidirectional_data <- bidirectional_data %>%
  group_by(query, subject) %>%
  summarise(ANI = max(ANI)) %>%  # Fill matrix with the higher score
  ungroup()
strains <- unique(c(bidirectional_data$query, bidirectional_data$subject))
ani_matrix <- matrix(NA, 
                     nrow = length(strains),
                     ncol = length(strains),
                     dimnames = list(strains, strains))
for (i in 1:nrow(bidirectional_data)) {
  q <- bidirectional_data$query[i]
  s <- bidirectional_data$subject[i]
  ani <- bidirectional_data$ANI[i]
  ani_matrix[q, s] <- ani
}
diag(ani_matrix) <- 100

# Function: determine typestrain for valid species
determine_prefix <- function(names) {
  prefixes <- character(length(names))
  for (i in seq_along(names)) {
    if (startsWith(names[i], "A2165.fna")) {
      prefixes[i] <- "F.duncaniae"
    } else if (startsWith(names[i], "GCA_036632915.2_ASM3663291v2_genomic.fna")) {
      prefixes[i] <- "F.taiwanense"
    } else if (startsWith(names[i], "APC922_41-1.fna")) {
      prefixes[i] <- "F.hattorii"
    } else if (startsWith(names[i], "CLAJMH7B.combined.fa")) {
      prefixes[i] <- "F.faecis"
    } else if (startsWith(names[i], "CLAAAH175.combined.fa")) {
      prefixes[i] <- "F.tardum"
    } else if (startsWith(names[i], "CLAAAH281.combined.fa")) {
      prefixes[i] <- "F.intestinale"
    } else if (startsWith(names[i], "CM04-06A.fna")) {
      prefixes[i] <- "F.longum"
    } else if (startsWith(names[i], "AF52-21.fna")) {
      prefixes[i] <- "F.butyricigenerans"
    } else if (startsWith(names[i], "ATCC_27768.fna")) {
      prefixes[i] <- "F.prausnitzii"
    } else if (startsWith(names[i], "GCA_023347535.1_ASM2334753v1_genomic.fna")) {
      prefixes[i] <- "F.wellingii"
    } else if (startsWith(names[i], "GCF_002549775.1_ASM254977v1_genomic.fna")) {
      prefixes[i] <- "F.langellae"
    } else {
      prefixes[i] <- "Other"
    }
  }
  return(prefixes)
}
# Color mapping for typestrain
prefix_colors <- list(
  Prefix = c(
    "F.longum" = '#FF7256',
    "F.prausnitzii" = '#FFA500',
    "F.duncaniae" = "#FDD80A",
    "F.faecis" = "#556B2F",
    "F.hattorii" = "#A1D99B",
    "F.butyricigenerans" = "#7FFFd4",
    "F.wellingii" = "#8EE5EE",
    "F.langellae" = "#6495ED",
    "F.taiwanense" = "#B23AEE",
    "F.tardum" = "#8B008B",
    "F.intestinale" = "#FF1493",
    "Other" = "#CCCCCC"
  )
)

# Extract row and column names
row_names <- rownames(ani_matrix)
col_names <- colnames(ani_matrix)
# Create row and column annotation
row_anno <- data.frame(Prefix = determine_prefix(row_names))
rownames(row_anno) <- row_names
col_anno <- data.frame(Prefix = determine_prefix(col_names))
rownames(col_anno) <- col_names
# Segmented color matching
break_points <- c(0, 95, 100)  # 0-95-100
yellow_gradient <- colorRampPalette(c("#BF5B17","#FDC086","#FFFF99"))(95)
blue_gradient <- colorRampPalette(c("#BEAED4","#386CB0"))(5)
color_map <- c(yellow_gradient, blue_gradient)
color_breaks <- c(seq(0, 95, length.out = 96), seq(95.1, 100, length.out = 5))

# Draw and save heatmap
library(pheatmap)
library(ggplot2)
p <- pheatmap(
  ani_matrix,
  annotation_row = row_anno,
  annotation_col = col_anno,
  annotation_colors = prefix_colors,
  color = color_map,
  breaks = color_breaks,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "average",
  main = "ANI Heatmap of Genus Faecalibacterium with 95 Cutoff",
  legend = TRUE,
  legend_breaks = c(47.5, 97.5),
  legend_labels = c("<95", "≥95"),
  legend_title = "ANI Value",
  border_color = ifelse(ani_matrix >= 95, "#001A4C", "#FFD700"),
  silent = TRUE
)
ggsave(
  filename = "ANIheatmap_Genus.pdf",
  plot = p$gtable,
  device = "pdf",
  width = 16,
  height = 12,
  units = "in",
  dpi = 600,
  useDingbats = FALSE
)