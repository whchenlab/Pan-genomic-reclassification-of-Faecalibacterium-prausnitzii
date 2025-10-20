# 导入及预处理
ani_data <- read.table("ani_matrix_genus.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)[,1:3]
colnames(ani_data) <- c("query", "subject", "ANI")
ani_data[, 1:2] <- lapply(ani_data[, 1:2], function(x) sub("\\.fna$", "", x))
bidirectional_data <- rbind(
  ani_data,
  data.frame(
    query = ani_data$subject,
    subject = ani_data$query,
    ANI = ani_data$ANI
  )
)


# 处理不一致的双向记录并填充
library(dplyr)
bidirectional_data <- bidirectional_data %>%
  group_by(query, subject) %>%
  summarise(ANI = max(ANI)) %>%  # 取高值
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

# Typestrain判断函数
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
      prefixes[i] <- "Other"  # 非Typestrain
    }
  }
  return(prefixes)
}

# 颜色映射
prefix_colors <- list(
  Prefix = c(
    "F.duncaniae" = "#FF5733",
    "F.taiwanense" = "#33FF57",
    "F.hattorii" = "#3357FF",
    "F.faecis" = "#FF33F6",
    "F.tardum" = "#33FFF6",
    "F.intestinale" = "#F6FF33",
    "F.longum" = "#FF3333",
    "F.butyricigenerans" = "#33FF33",
    "F.prausnitzii" = "#F633FF",
    "F.wellingii" = "#FFA500",
    "F.gallinarum" = "#8A2BE2",
    "Other" = "#CCCCCC"
  )
)

# 提取行名和列名
row_names <- rownames(ani_matrix)
col_names <- colnames(ani_matrix)
# 创建行和列注释数据框
row_anno <- data.frame(Prefix = determine_prefix(row_names))
rownames(row_anno) <- row_names
col_anno <- data.frame(Prefix = determine_prefix(col_names))
rownames(col_anno) <- col_names
# 分段配色
break_points <- c(0, 95, 100)  # 断点：0-95-100
yellow_gradient <- colorRampPalette(c("#BF5B17","#FDC086","#FFFF99"))(95)
blue_gradient <- colorRampPalette(c("#BEAED4","#386CB0"))(5)
color_map <- c(yellow_gradient, blue_gradient)
color_breaks <- c(seq(0, 95, length.out = 96), seq(95.1, 100, length.out = 5))

# 绘制热图
library(pheatmap)
library(ggplot2)
# 生成热图对象，并赋值给变量
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

# 保存为PDF
ggsave(
  filename = "ANIheatmap_Genus.pdf",
  plot = p$gtable,
  device = "pdf",
  width = 8,
  height = 6.4,
  units = "in",
  dpi = 600,
  useDingbats = FALSE
)

