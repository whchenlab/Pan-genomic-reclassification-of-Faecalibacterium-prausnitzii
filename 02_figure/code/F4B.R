library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(grid)

# 定义颜色方案
col.scheme.heatmap <- c('#F7FBFF', 'steelblue1', '#08306B')

# 函数：创建混淆矩阵并绘制热图
plot_confusion_matrix <- function(data, dataset_name) {
  # 转换为矩阵
  matrix_data <- matrix(
    data = data,
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Faecalibacterium", "Non-Faecalibacterium"),  # 真实类别
      c("Faecalibacterium", "Non-Faecalibacterium")   # 预测类别
    )
  )
  
  # 数据重塑
  df_counts <- melt(matrix_data, 
                    varnames = c("True Class", "Predicted Class"), 
                    value.name = "Count")
  
  # 确定文字颜色阈值
  max_count <- max(matrix_data)
  text_threshold <- ifelse(max_count == 0, 0, max_count / 2)
  
  # 绘制热图
  ggplot(df_counts, aes(x = `Predicted Class`, y = `True Class`)) +
    geom_tile(aes(fill = Count), color = "white", linewidth = 0.5) +
    geom_text(aes(label = Count), 
              color = ifelse(df_counts$Count > text_threshold, "white", "black"),
              size = 4, fontface = "bold") +
    scale_fill_gradientn(colours = col.scheme.heatmap, limits = c(0, 1740)) + 
    labs(
      title = dataset_name,
      x = "Predicted Class",
      y = "True Class",
      fill = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10, margin = margin(b = 5)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
      axis.text.y = element_text(size = 8, face = "italic"),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
    )
}

# 四个数据集的数据
ncbi_data <- c(482, 0, 0, 0)
hibc_data <- c(11, 0, 0, 329)
allfna2_data <- c(136, 0, 0, 0)
cgr2_data <- c(11, 0, 0, 1740)

# 生成四个混淆矩阵图
p1 <- plot_confusion_matrix(ncbi_data, "NCBI")
p2 <- plot_confusion_matrix(hibc_data, "HIBC")
p3 <- plot_confusion_matrix(allfna2_data, "BGI")
p4 <- plot_confusion_matrix(cgr2_data, "CGR2")

# 创建共享图例
legend <- get_legend(
  p1 + theme(legend.position = "bottom", 
             legend.key.width = unit(2, "cm"),
             legend.title = element_text(face = "bold", size = 10),
             legend.text = element_text(size = 9))
)

# 组合四个子图
subplots <- arrangeGrob(p1, p2, p3, p4, ncol = 2)

# 创建大标题文本
main_title <- textGrob(
  "Genus-Level Classification Performance for Faecalibacterium",
  gp = gpar(fontface = "bold", fontsize = 14)
)

# 组合所有元素（标题 + 子图 + 图例）
combined_plot <- grid.arrange(
  main_title,
  subplots,
  legend,
  nrow = 3,
  heights = c(0.5, 10, 1)
)

# 显示组合图
print(combined_plot)

