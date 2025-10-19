# 读取数据，更改文件路径
ani <- read.table("ani_matrix_genus.txt", sep = "\t")

# 加载所需库
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# 长表转宽表，处理列
ani_t <- ani[, -c(4, 5)]  # 移除不需要的列
ani_long <- ani_t %>%
  spread(key = V1, value = V3)  # 转换为宽表格式

# 设置行名并清理名称（去除.fna后缀）
row.names(ani_long) <- ani_long$V2
ani_long$V2 <- NULL
rownames(ani_long) <- gsub(".fna", "", rownames(ani_long))
colnames(ani_long) <- gsub(".fna", "", colnames(ani_long))

# 定义函数：根据菌株名前缀获取对应的typestrain名称
get_typestrain_name <- function(name) {
  if (startsWith(name, "A2165")) {
    return("F.duncaniae")
  } else if (startsWith(name, "GCA_036632915")) {
    return("F.taiwanense")
  } else if (startsWith(name, "APC922_41-1")) {
    return("F.hattorii")
  } else if (startsWith(name, "CLAJMH7B.combined")) {
    return("F.faecis")
  } else if (startsWith(name, "CLAAAH175.combined")) {
    return("F.tardum")
  } else if (startsWith(name, "CLAAAH281.combined")) {
    return("F.intestinale")
  } else if (startsWith(name, "scaffold_CM04-06A")) {
    return("F.longum")
  } else if (startsWith(name, "AF52-21")) {
    return("F.butyricigenerans")
  } else if (startsWith(name, "ATCC_27768")) {
    return("F.prausnitzii")
  } else if (startsWith(name, "F.wellingii_HTF-F")) {
    return("F.wellingii")
  } else if (startsWith(name, "F.gallinarum-JCM-17207")) {
    return("F.gallinarum")
  } else {
    return(NA)  # 非typestrain返回NA
  }
}

# 筛选出所有typestrain的名称
all_names <- union(rownames(ani_long), colnames(ani_long))
typestrain_flags <- !is.na(sapply(all_names, get_typestrain_name))
typestrain_names <- all_names[typestrain_flags]

# 提取仅包含typestrain的矩阵
ani_typestrain <- ani_long[
  rownames(ani_long) %in% typestrain_names,
  colnames(ani_long) %in% typestrain_names
]

# 转换为矩阵格式以便处理
ani_typestrain_mat <- as.matrix(ani_typestrain)

# 处理双相记录不一致：取对称位置的较高值
ani_typestrain_mat <- pmax(ani_typestrain_mat, t(ani_typestrain_mat))

# 替换行名和列名为对应的typestrain名称
name_mapping <- sapply(rownames(ani_typestrain_mat), get_typestrain_name)
rownames(ani_typestrain_mat) <- name_mapping
colnames(ani_typestrain_mat) <- name_mapping

# 准备热图参数
heatmap_colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)

# 处理显示数值：只显示大于95的值
ani_display <- round(ani_typestrain_mat, 1)
ani_display[ani_display < 95] <- NA
display_numbers_matrix <- ani_display
display_numbers_matrix[is.na(display_numbers_matrix)] <- ""

# 绘制并保存热图
pdf("ANIheatmap_Typestrain.pdf", width = 8, height = 6)
pheatmap(ani_typestrain_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         fontsize_row = 8, 
         fontsize_col = 8,
         color = heatmap_colors,
         main = "ANI Matrix between Faecalibacterium typestrains",
         display_numbers = display_numbers_matrix,
         number_color = "black",
         cellwidth = 20, 
         cellheight = 20,
         border_color = NA,
         treeheight_row = 30,
         treeheight_col = 30
)
dev.off()