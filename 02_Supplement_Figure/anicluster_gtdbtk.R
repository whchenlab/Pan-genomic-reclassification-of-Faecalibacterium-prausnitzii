# Load packages
library(ggplot2)
library(reshape2)

# Import data
data <- read.csv("gtdbtk_accuracy.csv", stringsAsFactors = FALSE)

# Define species
data$classification <- gsub("Faecalibacterium ", "F.", data$classification)
ani_order <- c(
  "F.longum", "F.prausnitzii", "F.duncaniae", "F.faecis",
  "F.un1", "F.un2", "F.hattorii", "F.un3", 
  "F.butyricigenerans", "F.un4", "F.wellingii", "F.langellae",
  "F.un5", "F.un6", "F.un7", "F.un8", "F.un9", "F.un10",
  "F.un11", "F.un12", "F.un13", "F.un14", "F.un15", "F.un16",
  "F.un17", "F.un18", "F.un19", "F.un20", "single"
)
ani_order <- intersect(ani_order, unique(data$ANIcluster))

# Create matrix: ANIcluster * classification
unique_class <- unique(data$classification)
count_matrix <- matrix(0, nrow = length(ani_order), ncol = length(unique_class),
                       dimnames = list(ani_order, unique_class))
for(i in 1:nrow(data)) {
  ani <- data$ANIcluster[i]
  class <- data$classification[i]
  count <- data$count[i]
  
  if(ani %in% ani_order) {
    count_matrix[ani, class] <- count
  }
}
col_nonzero <- colSums(count_matrix > 0)
col_positions <- apply(count_matrix > 0, 2, function(x) {
  if(any(x)) {
    which(x)[1]
  } else {
    Inf
  }
})

# Set order
col_total <- colSums(count_matrix)
col_info <- data.frame(
  class = colnames(count_matrix),
  first_pos = col_positions,
  nonzero_count = col_nonzero,
  total = col_total
)
col_info <- col_info[order(col_info$first_pos, -col_info$nonzero_count, -col_info$total), ]
class_order <- as.character(col_info$class)
count_matrix <- count_matrix[, class_order, drop = FALSE]
df_counts <- melt(count_matrix, 
                  varnames = c("ANIcluster", "Classification"), 
                  value.name = "Count")
df_counts <- df_counts[df_counts$Count > 0, ]
df_counts$ANIcluster <- factor(df_counts$ANIcluster, levels = ani_order)
df_counts$Classification <- factor(df_counts$Classification, levels = class_order)

# Set color mapping
col.scheme.heatmap <- c('#F7FBFF', 'steelblue1', '#08306B')

# Create plot
p1 <- ggplot(df_counts, aes(x = Classification, y = ANIcluster)) +
  geom_tile(aes(fill = Count), color = "white", linewidth = 0.5) +
  geom_text(aes(label = Count), 
            color = ifelse(df_counts$Count > max(df_counts$Count)/2, "white", "black"), 
            size = 2, fontface = "bold") +
  scale_fill_gradientn(colours = col.scheme.heatmap, 
                       limits = c(0, max(count_matrix)),
                       na.value = "white") +
  labs(
    title = " ",
    x = "GTDB-TK Classification",
    y = "ANI Cluster",
    fill = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
    axis.text.y = element_text(size = 8, face = "italic"),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(colour = "lightgray", fill = NA, linewidth = 1),
    plot.background = element_rect(colour = "lightgray", fill = NA, linewidth = 1)
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE)

# Print and save
print(p1)
ggsave("anicluster_gtdbtk.pdf", p1, width = 10, height = 6, dpi = 300)
