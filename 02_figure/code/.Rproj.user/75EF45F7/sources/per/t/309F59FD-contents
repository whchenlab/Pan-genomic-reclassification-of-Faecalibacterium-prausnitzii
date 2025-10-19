# 预处理：确保双向比对记录一致
ani_data <- read.table("ani_matrix_genus.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)[,1:3]
colnames(ani_data) <- c("query", "subject", "ANI")

# 去除.fna后缀
ani_data[, 1:2] <- lapply(ani_data[, 1:2], function(x) sub("\\.fna$", "", x))

# 生成双向一致的比对数据
bidirectional_data <- rbind(
  ani_data,
  data.frame(
    query = ani_data$subject,
    subject = ani_data$query,
    ANI = ani_data$ANI
  )
)

# 按(query, subject)分组，取最高ANI值（处理不一致的双向记录）
library(dplyr)
bidirectional_data <- bidirectional_data %>%
  group_by(query, subject) %>%
  summarise(ANI = max(ANI)) %>%  # 此处从mean改为max
  ungroup()

# 创建菌株列表并构建矩阵
strains <- unique(c(bidirectional_data$query, bidirectional_data$subject))
ani_matrix <- matrix(NA, 
                     nrow = length(strains),
                     ncol = length(strains),
                     dimnames = list(strains, strains))

# 填充矩阵（此时数据已双向一致）
for (i in 1:nrow(bidirectional_data)) {
  q <- bidirectional_data$query[i]
  s <- bidirectional_data$subject[i]
  ani <- bidirectional_data$ANI[i]
  ani_matrix[q, s] <- ani
}

# 设置对角线为100
diag(ani_matrix) <- 100

# 前缀判断函数（保持不变）
determine_prefix <- function(names) {
  prefixes <- character(length(names))
  for (i in seq_along(names)) {
    if (startsWith(names[i], "A2165") || 
        startsWith(names[i], "GCA_032971305") || 
        startsWith(names[i], "GCA_010509575")) {
      prefixes[i] <- "F.duncaniae"
    } else if (startsWith(names[i], "GCA_036632915")) {
      prefixes[i] <- "F.taiwanense"
    } else if (startsWith(names[i], "APC922_41-1")) {
      prefixes[i] <- "F.hattorii"
    } else if (startsWith(names[i], "CLAJMH7B.combined")) {
      prefixes[i] <- "F.faecis"
    } else if (startsWith(names[i], "CLAAAH175.combined")) {
      prefixes[i] <- "F.tardum"
    } else if (startsWith(names[i], "CLAAAH281.combined")) {
      prefixes[i] <- "F.intestinale"
    } else if (startsWith(names[i], "scaffold_CM04-06A")) {
      prefixes[i] <- "F.longum"
    } else if (startsWith(names[i], "AF52-21")) {
      prefixes[i] <- "F.butyricigenerans"
    } else if (startsWith(names[i], "ATCC_27768")) {
      prefixes[i] <- "F.prausnitzii"
    } else if (startsWith(names[i], "F.wellingii_HTF-F")) {
      prefixes[i] <- "F.wellingii"
    } else if (startsWith(names[i], "F.gallinarum-JCM-17207")) {
      prefixes[i] <- "F.gallinarum"
    } else {
      prefixes[i] <- "Other"  # 其他情况
    }
  }
  return(prefixes)
}

# 颜色映射（保持不变）
prefix_colors <- list(
  Prefix = c(
    "F.duncaniae" = "#FF5733",  # 橙色
    "F.taiwanense" = "#33FF57",  # 浅绿色
    "F.hattorii" = "#3357FF",  # 深蓝色
    "F.faecis" = "#FF33F6",  # 紫色
    "F.tardum" = "#33FFF6",  # 青色
    "F.intestinale" = "#F6FF33",  # 黄色
    "F.longum" = "#FF3333",  # 红色
    "F.butyricigenerans" = "#33FF33",  # 绿色
    "F.prausnitzii" = "#F633FF",  # 粉色
    "F.wellingii" = "#FFA500",   # 橙色
    "F.gallinarum" = "#8A2BE2",  # 靛蓝色
    "Other" = "#CCCCCC"  # 灰色
  )
)

# 提取行名和列名（保持不变）
row_names <- rownames(ani_matrix)
col_names <- colnames(ani_matrix)
# 创建行和列注释数据框
row_anno <- data.frame(Prefix = determine_prefix(row_names))
rownames(row_anno) <- row_names
col_anno <- data.frame(Prefix = determine_prefix(col_names))
rownames(col_anno) <- col_names

# 分段的黄蓝配色（保持不变）
break_points <- c(0, 95, 100)  # 断点：0-95-100
yellow_gradient <- colorRampPalette(c("#BF5B17","#FDC086","#FFFF99"))(95)  # 0-95的黄色渐变
blue_gradient <- colorRampPalette(c("#BEAED4","#386CB0"))(5)    # 95-100的蓝色渐变
color_map <- c(yellow_gradient, blue_gradient)     # 合并两个颜色渐变
color_breaks <- c(seq(0, 95, length.out = 96), seq(95.1, 100, length.out = 5))     # 生成对应的断点

# 绘制热图
library(pheatmap)
library(ggplot2)
# 先生成热图对象并赋值给变量
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
  silent = TRUE  # 静默模式，不在控制台输出图形
)

# 使用ggsave保存为PDF，指定宽高比例
ggsave(
  filename = "ANIheatmap_Genus.pdf",  # 保存的文件名
  plot = p$gtable,               # pheatmap返回的gtable对象
  device = "pdf",                # 设备类型
  width = 8,                    # 宽度（英寸）
  height = 6.4,                    # 高度（英寸），10:8的比例
  units = "in",                  # 单位
  dpi = 600,                     # 分辨率
  useDingbats = FALSE            # 避免特殊符号问题
)

