library(ggplot2)
library(dplyr)

# 1. 数据预处理
strain_cols <- 15:ncol(matrix)
presence_matrix <- matrix[, strain_cols, drop = FALSE]
# 转换为存在/缺失矩阵
presence_matrix[presence_matrix != "" & !is.na(presence_matrix)] <- 1
presence_matrix[presence_matrix == "" | is.na(presence_matrix)] <- 0
presence_matrix <- as.data.frame(lapply(presence_matrix, as.numeric))
# 计算基因出现次数和总菌株数
gene_counts <- rowSums(presence_matrix)
n_genomes <- ncol(presence_matrix)

# 2. 基因类型统计（固定顺序，包含基因数量）
pangenome_stats <- data.frame(
  category = factor(
    c("core", "softcore", "shell", "cloud"),
    levels = c("core", "softcore", "shell", "cloud")
  ),
  count = c(
    sum(gene_counts >= 0.95 * n_genomes),
    sum(gene_counts >= 0.94 * n_genomes & gene_counts < 0.95 * n_genomes),
    sum(gene_counts >= 0.15 * n_genomes & gene_counts < 0.94 * n_genomes),
    sum(gene_counts < 0.15 * n_genomes)
  ),
  strain_range = c(
    paste0("≥", round(0.95 * n_genomes), " strains"),
    paste0(round(0.94 * n_genomes), "–", round(0.95 * n_genomes), " strains"),
    paste0(round(0.15 * n_genomes), "–", round(0.94 * n_genomes), " strains"),
    paste0("<", round(0.15 * n_genomes), " strains")
  )
) %>%
  # 组合图例文本（基因类型 + 菌株范围 + 基因数量）
  mutate(
    legend_text = paste0(
      category, "\n",          # 基因类型
      strain_range, "\n",      # 菌株范围
      "n=", count              # 基因数量
    )
  )

# 3. 饼图参数
pie_x <- 0.5    # 饼图中心x坐标（0-1）
pie_y <- 0.55   # 饼图中心y坐标（0-1）
pie_r <- 0.5    # 饼图半径

# 4. 定义高区分度颜色方案
pie_colors <- c(
  "#e41a1c",  # 红色（core）
  "#4daf4a",  # 绿色（softcore）
  "#377eb8",  # 蓝色（shell）
  "#984ea3"   # 紫色（cloud）
)

# 5. 绘制频率直方图（主图）
p <- ggplot(data.frame(gene_counts), aes(x = gene_counts)) +
  geom_histogram(
    binwidth = 1,
    fill = "#f0f0f0",
    color = "black",
    alpha = 0.7
  ) +
  # 图例显示完整信息（类型+范围+数量）
  scale_fill_manual(
    values = pie_colors,
    name = "Gene Categories",
    labels = pangenome_stats$legend_text
  ) +
  labs(
    x = "Number of genomes with gene",
    y = "Number of genes",
    title = "Pangenome Distribution with Pie Chart"
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

# 6. 计算饼图坐标
x_lim <- layer_scales(p)$x$get_limits()
y_lim <- layer_scales(p)$y$get_limits()
x_length <- x_lim[2] - x_lim[1]
y_length <- y_lim[2] - y_lim[1]
aspect_ratio <- y_length / x_length

pie_x_actual <- x_lim[1] + pie_x * x_length
pie_y_actual <- y_lim[1] + pie_y * y_length
pie_r_actual <- pie_r * x_length

# 7. 创建饼图
pie_plot <- ggplot(pangenome_stats, aes(x = 0, y = count, fill = category)) +
  geom_col(width = 1, color = "white", size = 0.7) +
  coord_polar("y", start = 0) +
  xlim(-1, 1) +
  ylim(0, sum(pangenome_stats$count)) +
  scale_fill_manual(values = pie_colors, guide = "none") +
  theme_void()

# 8. 嵌入饼图到主图
pie_grob <- ggplotGrob(pie_plot)
p <- p + annotation_custom(
  grob = pie_grob,
  xmin = pie_x_actual - pie_r_actual,
  xmax = pie_x_actual + pie_r_actual,
  ymin = pie_y_actual - pie_r_actual * aspect_ratio,
  ymax = pie_y_actual + pie_r_actual * aspect_ratio
)

# 9. 图例显示
p <- p + geom_point(
  data = pangenome_stats,
  aes(x = -Inf, y = -Inf, fill = category),
  shape = 21,
  size = 6,
  color = "black",
  show.legend = TRUE
)

