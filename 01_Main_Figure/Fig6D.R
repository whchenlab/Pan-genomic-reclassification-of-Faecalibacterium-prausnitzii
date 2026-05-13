# Load required packages
library(ggplot2)
library(reshape2)
library(dplyr)

# Import and convert data
confusion_matrix_binary <- matrix(
  data = c(
    8913, 0,
    0, 2069
  ),
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("Faecalibacterium", "Non-Faecalibacterium"),  #True
    c("Faecalibacterium", "Non-Faecalibacterium")  #Predicted
  )
)
df_counts_binary <- melt(confusion_matrix_binary, 
                         varnames = c("True Class", "Predicted Class"), 
                         value.name = "Count")

# Create heatmap
col.scheme.heatmap <- c('#F7FBFF', 'steelblue1', '#08306B')
p_binary <- ggplot(df_counts_binary, aes(x = `Predicted Class`, y = `True Class`)) +
  geom_tile(aes(fill = Count), color = "white", linewidth = 0.5) +
  geom_text(aes(label = Count), 
            color = ifelse(df_counts_binary$Count > 50, "white", "black"), 
            size = 5, fontface = "bold") +
  scale_fill_gradientn(colours = col.scheme.heatmap, limits = c(0, max(confusion_matrix_binary))) +
  labs(title = " ", x = "Predicted Class", y = "True Class", fill = "Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 10, face = "italic"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# Save as pdf
pdf(
  file = "confusion_matrxi_genus_70.pdf",
  width = 6,
  height = 5,
)
print(p_binary)
dev.off()
