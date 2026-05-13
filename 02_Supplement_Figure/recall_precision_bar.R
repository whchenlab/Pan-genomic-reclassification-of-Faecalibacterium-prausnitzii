# Load required packages
library(ggplot2)
library(reshape2)
library(scales)
library(patchwork)
library(dplyr)
library(tidyr)
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
# 5. Calculate recall,precision and accuracy
# --------------------------
# Recall
recall_values <- diag(confusion_matrix) / rowSums(confusion_matrix)
recall_percent <- recall_values * 100

# Precision
precision_values <- diag(confusion_matrix) / colSums(confusion_matrix)
precision_percent <- precision_values * 100

# Overall accuracy
correct_predictions <- sum(diag(confusion_matrix))
total_predictions <- sum(confusion_matrix)
overall_accuracy <- correct_predictions / total_predictions
overall_accuracy_percent <- overall_accuracy * 100

# Create dataframe of recall,precision and accuracy
plot_data <- data.frame(
  Classifier = factor(rep(names(recall_percent), 2), levels = class_order),
  Metric = factor(c(rep("Recall", length(recall_percent)), 
                    rep("Precision", length(precision_percent))),
                  levels = c("Recall", "Precision")),
  Value = c(recall_percent, precision_percent)
)
plot_data$Value <- pmin(pmax(plot_data$Value, 0), 100)

# --------------------------
# 6. Draw bar plot
# --------------------------
color <- c("#006D2C","#74C476")

p_accuracy <- ggplot(plot_data, aes(x = Classifier, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.75),
           width = 0.7,
           color = "black",
           alpha = 0.9) +
  # Add line to show overall accuracy
  geom_hline(yintercept = overall_accuracy_percent, 
             color = "darkred", 
             linetype = "dashed", 
             linewidth = 1.2,
             alpha = 0.8) +
  annotate("text", 
           x = length(class_order)/2 + 0.5,  # 中间位置
           y = overall_accuracy_percent + max(plot_data$Value) * 0.03,  # 略高于虚线
           label = paste0("Overall Accuracy: ", round(overall_accuracy_percent, 1), "%"),
           color = "darkred",
           size = 4,
           fontface = "bold",
           hjust = 0.5) +
  # Add value on bars
  geom_text(aes(label = paste0(round(Value, 1), "%"), 
                group = Metric),
            position = position_dodge(width = 0.75),
            vjust = 2,
            size = 1.5,
            fontface = "bold",
            color = "white") +
  # Set color for bars
  scale_fill_manual(
    name = "Group",
    values = c("Recall" = color[1],  
               "Precision" = color[2]),  
    labels = c("Recall", "Precision")
  ) +
  scale_y_continuous(
    name = "value",
    limits = c(0, max(plot_data$Value, overall_accuracy_percent) * 1.05),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(
    name = "Species"
  ) +
  labs(
    title = " "
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

print(p_accuracy)
ggsave("recall_precision_bar_90.pdf", p_accuracy, width = 6, height = 5, dpi = 300)
