library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(patchwork)

# --------------------------
# 1. 数据处理
# --------------------------
summary$Actual_short <- NA
# 处理F.prausnitzii类别（包含Fp_、F.tardum、F.intestinale前缀）
fp_pattern <- "^(Fp_|F.tardum|F.intestinale|F.prausnitzii)"
summary$Actual_short[grepl(fp_pattern, summary$Genome)] <- "Fp"
# 处理F.duncaniae类别（包含Fd_、F.duncaniae前缀）
fd_pattern <- "^(Fd_|F.duncaniae)"
summary$Actual_short[grepl(fd_pattern, summary$Genome) & is.na(summary$Actual_short)] <- "Fd"
# 处理F.longum类别（Fl_前缀）
fl_pattern <- "^Fl_"
summary$Actual_short[grepl(fl_pattern, summary$Genome) & is.na(summary$Actual_short)] <- "Fl"
# 剩余的都归为Other
summary$Actual_short[is.na(summary$Actual_short)] <- "Other"

# --------------------------
# 2. 构建混淆矩阵 - 四分类版本
# --------------------------
name_mapping <- c(
  "Fp" = "F.prausnitzii",
  "Fd" = "F.duncaniae",
  "Fl" = "F.longum",
  "Other" = "F.other"
)

full_name_order <- c("F.prausnitzii", "F.duncaniae", "F.longum", "F.other")
class_order_short <- c("Fp", "Fd", "Fl", "Other")
summary$Predicted_short <- factor(summary$Type, levels = class_order_short)

confusion_matrix <- table(
  `True Class` = factor(summary$Actual_short, levels = class_order_short),
  `Predicted Class` = summary$Predicted_short
)

# --------------------------
# 3. 绘制混淆矩阵热图
# --------------------------
df_counts <- melt(confusion_matrix, 
                  varnames = c("True Class", "Predicted Class"), 
                  value.name = "Count")

df_counts$`True Class` <- factor(
  name_mapping[as.character(df_counts$`True Class`)],
  levels = full_name_order
)
df_counts$`Predicted Class` <- factor(
  name_mapping[as.character(df_counts$`Predicted Class`)],
  levels = full_name_order
)

col.scheme.heatmap <- c('#F7FBFF', 'steelblue1', '#08306B')

# 热图
p1 <- ggplot(df_counts, aes(x = `Predicted Class`, y = `True Class`)) +
  geom_tile(aes(fill = Count), color = "white", linewidth = 0.5) +
  geom_text(aes(label = Count), 
            color = ifelse(df_counts$Count > 50, "white", "black"), 
            size = 4, fontface = "bold") +
  scale_fill_gradientn(colours = col.scheme.heatmap, 
                       limits = c(0, max(confusion_matrix))) +
  labs(
    title = "Confusion Matrix",
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
    panel.grid = element_blank()
  )

# --------------------------
# 4. 计算并绘制Precision和Recall
# --------------------------
# 总体准确率
overall_accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

# 每个类别的Precision和Recall
metrics <- data.frame(
  Class = full_name_order,
  Precision = diag(confusion_matrix) / colSums(confusion_matrix),
  Recall = diag(confusion_matrix) / rowSums(confusion_matrix)
)

metrics[is.na(metrics$Precision), "Precision"] <- 0
metrics[is.na(metrics$Recall), "Recall"] <- 0
metrics_melt <- melt(metrics, id.vars = "Class", variable.name = "Metric", value.name = "Value")
metrics_melt$Class <- factor(metrics_melt$Class, levels = full_name_order)

# 条形图
p2 <- ggplot(metrics_melt, aes(x = Class, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = overall_accuracy, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = length(full_name_order)/2, y = overall_accuracy + 0.02, 
           label = paste("Overall Accuracy =", round(overall_accuracy, 4)), 
           color = "red", fontface = "bold", size = 4) +
  scale_y_continuous(limits = c(0, 1.05), oob = rescale_none) +
  scale_fill_manual(values = c("#3274A1", "#E1812C"), 
                    labels = c("Precision", "Recall")) +
  labs(
    title = "Precision and Recall",
    x = "Species",
    y = "Score",
    fill = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "italic"),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

