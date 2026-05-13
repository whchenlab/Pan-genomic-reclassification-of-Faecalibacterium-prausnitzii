# Load required packages
library(ggplot2)
library(reshape2)
library(scales)
library(patchwork)
summary_with_ani <- summary_with_ani_90

# --------------------------
# 1. Add true classified type by ANIcluster
# --------------------------

summary_with_ani$Actual <- NA
for (i in 1:nrow(summary_with_ani)) {
  ani_cluster <- as.character(summary_with_ani$ANIcluster[i])
  
  # ANIcluster to classified type
  actual_class <- switch(
    ani_cluster,
    "F.longum" = "F.longum",
    "F.prausnitzii" = "F.prausnitzii",
    "F.un16" = "F.prausnitzii",
    "F.un9" = "F.prausnitzii", 
    "F.un10" = "F.prausnitzii",
    "F.un5" = "F.prausnitzii",
    "F.duncaniae" = "F.duncaniae",
    "F.faecis" = "F.faecis",
    "F.un1" = "F.un1",
    "F.un2" = "F.un2",
    "F.hattorii" = "F.hattorii",
    "F.un3" = "F.un3",
    "F.butyricigenerans" = "F.butyricigenerans",
    "F.un4" = "F.un4",
    "F.wellingii" = "F.wellingii",
    "F.langellae" = "F.langellae",
    "Other"
  )
  
  # Other=F.else
  if (is.null(actual_class) || is.na(actual_class)) {
    if (grepl("^F\\.un", ani_cluster)) {
      summary_with_ani$Actual[i] <- "F.else"
    } else {
      summary_with_ani$Actual[i] <- "F.else"
    }
  } else if (actual_class == "Other") {
    summary_with_ani$Actual[i] <- "F.else"
  } else {
    summary_with_ani$Actual[i] <- actual_class
  }
}

# --------------------------
# 2. Add predicted classified type
# --------------------------

summary_with_ani$Predicted <- NA
# Threshold for F.else
for (i in 1:nrow(summary_with_ani)) {
  if (summary_with_ani$Mean_MScore[i] < 139) {
    summary_with_ani$Predicted[i] <- "F.else"
  } else {
    # Predicted type to classified type
    max_type <- as.character(summary_with_ani$Max_Type[i])
    predicted_class <- switch(
      max_type,
      "Fl" = "F.longum",
      "Fp" = "F.prausnitzii", 
      "Fd" = "F.duncaniae",
      "Ff" = "F.faecis",
      "Fua" = "F.un1",
      "Fub" = "F.un2",
      "Fh" = "F.hattorii",
      "Fuc" = "F.un3",
      "Fb" = "F.butyricigenerans",
      "Fud" = "F.un4",
      "Fw" = "F.wellingii",
      "Fla" = "F.langellae",
      "F.else"
    )

    if (is.null(predicted_class) || is.na(predicted_class)) {
      summary_with_ani$Predicted[i] <- "F.else"
    } else if (predicted_class == "F.else") {
      summary_with_ani$Predicted[i] <- "F.else"
    } else {
      summary_with_ani$Predicted[i] <- predicted_class
    }
  }
}

# --------------------------
# 3. Define classified type
# --------------------------
all_predicted <- unique(summary_with_ani$Predicted)
all_actual <- unique(summary_with_ani$Actual)
all_classes <- unique(c(all_predicted, all_actual))

# Set order
full_class_order <- c(
  "F.longum", "F.prausnitzii", "F.duncaniae", "F.faecis",
  "F.un1", "F.un2", "F.hattorii", "F.un3", 
  "F.butyricigenerans", "F.un4", "F.wellingii", "F.langellae",
  "F.else"
)

class_order <- intersect(full_class_order, all_classes)
summary_with_ani$Predicted <- factor(summary_with_ani$Predicted, levels = class_order)
summary_with_ani$Actual <- factor(summary_with_ani$Actual, levels = class_order)

# --------------------------
# 4. Create confusion matrix
# --------------------------
confusion_matrix <- table(
  `True Class` = summary_with_ani$Actual,
  `Predicted Class` = summary_with_ani$Predicted
)

# --------------------------
# 5. Calculate accuracy
# --------------------------
# Overall accuracy
correct_predictions <- sum(diag(confusion_matrix))
total_predictions <- sum(confusion_matrix)
accuracy <- correct_predictions / total_predictions
print(paste("\nOverall accuracy:", round(accuracy * 100, 2), "%"))

# Classified type accuracy
class_accuracy <- data.frame(
  TrueClass = rownames(confusion_matrix),
  Accuracy = diag(confusion_matrix) / rowSums(confusion_matrix),
  Total = rowSums(confusion_matrix)
)
class_accuracy$Accuracy_Label <- paste0(
  round(class_accuracy$Accuracy * 100, 1), "%"
)

accuracy_labels <- setNames(
  class_accuracy$Accuracy_Label,
  class_accuracy$TrueClass
)

# --------------------------
# 6. Draw heatmap
# --------------------------
# Convert matrix to dataframe
df_counts <- melt(confusion_matrix, 
                  varnames = c("True Class", "Predicted Class"), 
                  value.name = "Count")
df_counts_nonzero <- df_counts[df_counts$Count > 0, ]

# Create accuracy labels
accuracy_df <- data.frame(
  `True Class` = factor(names(accuracy_labels), levels = class_order),
  Accuracy = accuracy_labels,
  stringsAsFactors = FALSE
)
names(accuracy_df)[1] <- "True Class"

# Set color pattern
col.scheme.heatmap <- c('#F7FBFF', 'steelblue1', '#08306B')

# Draw heatmap
n_classes <- length(class_order)
expansion_factor <- 0.15
p1 <- ggplot(df_counts, aes(x = `Predicted Class`, y = `True Class`)) +
  geom_tile(aes(fill = Count), color = "white", linewidth = 0.5) +
  # Only show >0
  geom_text(data = df_counts_nonzero,
            aes(label = Count), 
            color = ifelse(df_counts_nonzero$Count > max(df_counts_nonzero$Count)/2, "white", "black"), 
            size = 4, fontface = "bold") +
  # Add accuracy labels
  geom_text(data = accuracy_df,
            aes(x = length(class_order) + 0.7,
                y = `True Class`,
                label = Accuracy),
            color = "darkred",
            size = 2.5,
            fontface = "bold",
            hjust = 0) +
  # Add title
  annotate("text",
           x = length(class_order) + 0.7,
           y = length(class_order) + 0.5,
           label = "",
           color = "darkred",
           size = 4,
           fontface = "bold",
           hjust = 0) +
  scale_fill_gradientn(colours = col.scheme.heatmap, 
                       limits = c(0, max(confusion_matrix)),
                       na.value = "white") +
  labs(
    title = "Confusion Matrix(90 Completeness Genomes)",
    x = "Predicted Class",
    y = "True Class",
    fill = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
    axis.text.y = element_text(size = 9, face = "italic"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(10, 80, 10, 10),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.ticks = element_line()
  ) +
  scale_x_discrete(drop = FALSE, expand = expansion(mult = c(0, expansion_factor))) +
  scale_y_discrete(drop = FALSE) +
  geom_vline(xintercept = length(class_order) + 0.35, 
             color = "gray70", linetype = "dashed", size = 0.5)
print(p1)
ggsave("confusion_matrix_90_90threshold.pdf", p1, width = 8, height = 6, dpi = 300)

