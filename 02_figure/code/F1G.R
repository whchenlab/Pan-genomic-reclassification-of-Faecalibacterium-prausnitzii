# 加载必要的包
library(stringr)
library(VennDiagram)
library(grid)

# 假设数据框为 eggnog_FP、eggnog_FD、eggnog_FL，均包含 eggNOG_OGs 列

# 1. 提取eggNOG_OGs中第一个@之前的内容
extract_og_id <- function(og_string) {
  ifelse(is.na(og_string), NA, 
         str_split_fixed(og_string, "@", n=2)[,1])
}

# 处理三个数据框并去重
fp_ogs <- unique(extract_og_id(eggnog_FP$eggNOG_OGs))
fd_ogs <- unique(extract_og_id(eggnog_FD$eggNOG_OGs))
fl_ogs <- unique(extract_og_id(eggnog_FL$eggNOG_OGs))

# 移除NA值
fp_ogs <- fp_ogs[!is.na(fp_ogs)]
fd_ogs <- fd_ogs[!is.na(fd_ogs)]
fl_ogs <- fl_ogs[!is.na(fl_ogs)]

# 2. 定义集合（用于韦恩图）
set_fp <- fp_ogs
set_fd <- fd_ogs
set_fl <- fl_ogs

# 3. 绘制韦恩图
venn.plot <- venn.diagram(
  x = list(
    "F. prausnitzii" = set_fp,
    "F. duncaniae" = set_fd,
    "F. longum" = set_fl
  ),
  filename = NULL,
  col = "transparent",
  fill = c("skyblue", "lightgreen", "lightpink"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.2,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("skyblue", "lightgreen", "lightpink"),
  cat.cex = 1.2,
  cat.fontfamily = "serif",
  cat.pos = 0,
  cat.dist = 0.02,
  rotation.degree = 270,
  margin = 0.1,
  main = "Overlap of Core Genes(95%) with eggNOG Annotation",
  main.cex = 1.8,
  main.fontfamily = "serif",
  main.fontface = "bold",
  sub = "Comparison of  F. prausnitzii, F. duncaniae and F. longum",
  sub.cex = 1.2,
  sub.fontfamily = "serif",
  sub.fontface = "italic"
)

# 4. 显示韦恩图
grid.newpage()
grid.draw(venn.plot)

# 5. 输出统计结果
cat("基于eggNOG_OGs的基因集合统计结果：\n")
cat(paste("三者共有的基因数量：", length(intersect(intersect(fp_ogs, fd_ogs), fl_ogs)), "\n", sep = ""))
cat(paste("FP和FD共有的基因数量（不含三者共有）：", length(setdiff(intersect(fp_ogs, fd_ogs), fl_ogs)), "\n", sep = ""))
cat(paste("FP和FL共有的基因数量（不含三者共有）：", length(setdiff(intersect(fp_ogs, fl_ogs), fd_ogs)), "\n", sep = ""))
cat(paste("FD和FL共有的基因数量（不含三者共有）：", length(setdiff(intersect(fd_ogs, fl_ogs), fp_ogs)), "\n", sep = ""))
cat(paste("FP特有的基因数量：", length(setdiff(fp_ogs, union(fd_ogs, fl_ogs))), "\n", sep = ""))
cat(paste("FD特有的基因数量：", length(setdiff(fd_ogs, union(fp_ogs, fl_ogs))), "\n", sep = ""))
cat(paste("FL特有的基因数量：", length(setdiff(fl_ogs, union(fp_ogs, fd_ogs))), "\n", sep = ""))
