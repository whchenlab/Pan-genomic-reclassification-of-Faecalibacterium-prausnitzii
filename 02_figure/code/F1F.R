library(ape)
library(ggtree)
library(ggplot2)

# 绘制进化树
# 参数设置
tree_file <- "mytree.newick"  # 进化树
output_file <- "phylogenetic_tree_with_labels.pdf"
label_size <- 3
tree_width <- 10
tree_height <- 8
# 读取绘图
tree <- read.tree(tree_file)
tree_plot <- ggtree(tree) +
  theme_tree() +
  # 添加菌株名
  geom_tiplab(size = label_size, align = TRUE, linesize = 0.2, color = "black") +
  labs(title = paste("进化树 (", length(tree$tip.label), " 个菌株)")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10, "mm")
  )
print(tree_plot)
# 保存为PDF
ggsave(
  filename = output_file,
  plot = tree_plot,
  device = "pdf",
  width = tree_width,
  height = tree_height,
  dpi = 300,
  limitsize = FALSE  # 允许图形尺寸超过默认限制
)

library(ape)
library(ggplot2)
library(ggtree)
library(tidyr)
library(dplyr)
library(cowplot)
library(reshape2)

# 定义268个菌株
strain_order_268 <- c(
  "GCA_022136525", "GCA_019061315", "GCA_022137305", "BIOML-B1", "GCA_048490295", 
  "GCA_048396715", "APC918_95b", "GCA_934825545", "GCA_934841625", "GCA_934882155", 
  "GCA_934835215", "GCA_965148815", "GCA_003434015", "AF29-11BH", "scaffold_AF29-7BH-C", 
  "GCA_003433905", "AF31-14AC", "GCA_048251395", "GCA_041324535", "GCA_022750585", 
  "GCA_022750565", "Fp1043", "Fp137", "GCA_902373685", "CLAJMH7.combined", "SL3_3", 
  "GCA_048458355", "GCA_048346755", "GCA_959027385", "GCA_048285995", "GCA_048363475", 
  "APC924_119", "GCA_004167405", "BIOML-B3", "BIOML-B18", "BIOML-B2", "BIOML-B17", 
  "BIOML-B4", "CLAAAH222.combined", "CNCMI4573", "GCA_048485695", "GCA_049532935", 
  "scaffold_AM54-7NS", "GCA_008680095", "GCA_003467805", "AM36-18BH", "GCA_047946905", 
  "GCA_041333395", "scaffold_AF32-9BH", "GCA_905200815", "GCA_041323705", "GCA_003465525", 
  "AF10-13", "GCA_021589125", "scaffold_AM100_B16A", "GCA_003434135", "AM37-13AC", 
  "GCA_937924275", "GCA_024463025", "GCA_024462975", "GCA_019061285", "GCA_048399835", 
  "GCA_958423435", "GCA_021626805", "GCA_048349215", "GCA_048205145", "GCA_024460915", 
  "GCA_048371115", "ATCC_27768", "ATCC_27766", "GCA_046884385", "BIOML-B9", "BIOML-B7", 
  "BIOML-B12", "BIOML-B11", "BIOML-B6", "BIOML-B5", "BIOML-B8", "BIOML-B14", "BIOML-B10", 
  "BIOML-B13", "GCA_958448145", "GCA_958445435", "GCA_041328415", "GCA_959600295", 
  "GCA_048404765", "GCA_048447705", "GCA_018784645", "GCA_041334015", "GCA_001406255", 
  "CNCMI4644", "CLAAAH175.combined", "OM04-11BH", "GCA_003482185", "P9238", "P9094", 
  "GCA_017377545", "CNCMI4574", "CNCMI4543", "GCA_032971445", "GCA_902362495", 
  "GCA_010509575", "A2165", "GCA_032971305", "GCA_036312825", "GCA_036301245", 
  "GCA_900066125", "GCA_036297725", "GCA_025993415", "GCA_048429305", "GCA_958452975", 
  "GCA_037344775", "GCA_020686905", "CLAAAH283.combined", "GCA_018784555", "GCA_018380685", 
  "GCA_032971385", "GCA_048292045", "P9225", "GCA_017377595", "GCA_965149985", 
  "GCA_032971355", "GCA_020687425", "CLAAAH223.combined", "GCA_018784565", "P9311", 
  "P9123", "CLAAAH16.combined", "Indica", "P9224", "GCA_958447025", "GCA_019415355", 
  "GCA_934882665", "GCA_934836715", "GCA_046995315", "BIOML-A1", "GCA_018375935", 
  "GCA_959608545", "GCA_937996555", "GCA_041404485", "GCA_020687145", "GCA_022750605", 
  "GCA_041405325", "GCA_019415165", "GCA_003434165", "AM39-7BH", "GCA_041404105", 
  "GCA_003434125", "AM42-11AC", "GCA_905204145", "GCA_049545105", "GCA_958477415", 
  "GCA_002550045", "GCA_041404045", "GCA_041405505", "Fp28", "Fp1", "Fp14", "GCA_003434175", 
  "AM33-14AC", "GCA_018379485", "GCA_041403865", "GCA_958352925", "GCA_014287915", 
  "P9391", "GCA_003433865", "AF36-11AT", "GCA_041405715", "GCA_048530285", "GCA_017377535", 
  "P9239", "P9241", "P9313", "P9240", "P9122", "GCA_048171345", "GCA_020687245", 
  "CLAAAH243.combined", "GCA_048339445", "GCA_959018785", "GCA_958441855", "GCA_958365925", 
  "GCA_018364195", "GCA_041331995", "GCA_041323105", "GCA_003293635", "APC942_18-1", 
  "GCA_049545765", "GCA_048511395", "GCA_017377665", "GCA_048286095", "GCA_965142725", 
  "GCA_048168425", "GCA_024463015", "GCA_019041425", "GCA_037459115", "GCA_003433995", 
  "AF32-8AC", "GCA_041326375", "GCA_048508875", "Fp360", "GCA_958434715", "GCA_958353575", 
  "GCA_041326645", "GCA_041328815", "GCA_041324655", "GCA_020687265", "CLAAAH236.combined", 
  "GCA_018378015", "GCA_905208395", "GCA_048199535", "GCA_025757765", "GCA_019415265", 
  "GCA_048302425", "GCA_048373655", "GCA_030839505", "GCA_959024715", "GCA_048475415", 
  "APC942_32-1", "GCA_048465695", "GCA_048218205", "GCA_047017165", "GCA_048463715", 
  "GCA_958421125", "GCA_030843165", "GCA_018368925", "GCA_018382675", "GCA_048486235", 
  "GCA_048376695", "Marseille-P9312", "GCA_018377455", "GCA_037322895", "GCA_047232625", 
  "GCA_047229995", "GCA_959028065", "GCA_048468875", "GCA_046982205", "GCA_046979045", 
  "GCA_048470115", "GCA_048349455", "scaffold_CM04-06A", "GCA_047232465", "GCA_018370385", 
  "GCA_019415285", "scaffold_AFF13-7A", "GCA_959601315", "GCA_048306095", "GCA_022740055", 
  "GCA_018381965", "FPSSTS7063_SV_a2", "JG_BgPS064", "GCA_022735375", "GCA_958415745", 
  "GCA_022715735", "GCA_047241425", "GCA_047238525", "GCA_047244245", "GCA_046972245", 
  "GCA_022721015", "GCA_048423265", "GCA_048397125", "GCA_018367245", "GCA_047249675", 
  "GCA_046948715", "L2_6", "GCA_046996585", "GCA_046977165", "GCA_048394125"
)
# 定义颜色分组边界
index1 <- which(strain_order_268 == "GCA_003482185")  # 第94个
index2_start <- which(strain_order_268 == "P9238")    # 第95个
index2_end <- which(strain_order_268 == "GCA_048171345")  # 第176个
index3 <- which(strain_order_268 == "GCA_020687245")  # 第177个
# 定义菌种-颜色映射
species_colors <- data.frame(
  Species = c("F.prausnitzii", "F.duncaniae", "F.longum"),
  Color = c("#984EA3", "#FF7F00", "#A65628"),
  stringsAsFactors = FALSE
)

# 参数设置
min_strains <- 78  # 最小菌株数
tree_file <- "mytree.newick"
output_format <- "pdf"
show_labels <- FALSE
# 读取进化树和基因存在缺失矩阵
tree <- read.tree(tree_file)
gene_data <- read.csv("gene_presence_absence.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 处理使两图顺序对应
# 1. 提取菌株列（第15列开始）
strain_columns <- gene_data[, 15:ncol(gene_data), drop = FALSE]
# 2. 0=缺失, 1=存在
binary_matrix <- apply(strain_columns, c(1,2), function(x) ifelse(nchar(x) > 0, 1, 0))
rownames(binary_matrix) <- gene_data$Gene
# 3. 筛选基因
gene_counts <- rowSums(binary_matrix)
filtered_matrix <- binary_matrix[gene_counts >= min_strains, , drop = FALSE]
# 检查是否有符合条件的基因
if (nrow(filtered_matrix) == 0) {
  stop(paste("未检测到存在于≥", min_strains, "个菌株中的基因，请调整参数或检查输入文件。"))
}
# 4. 按存在菌株数降序排序
gene_order <- order(rowSums(filtered_matrix), decreasing = TRUE)
filtered_matrix <- filtered_matrix[gene_order, , drop = FALSE]
# 5. 按指定的268个菌株顺序调整矩阵
common_strains <- intersect(colnames(filtered_matrix), strain_order_268)
if (length(common_strains) == 0) {
  stop("矩阵与指定的菌株列表中没有共同的菌株名，请检查名称是否匹配。")
}
# 6. 按指定顺序排列矩阵
final_order <- strain_order_268[strain_order_268 %in% common_strains]
filtered_matrix <- filtered_matrix[, final_order, drop = FALSE]
# 7. 同步调整进化树
tree <- keep.tip(tree, common_strains)
constraint_data <- data.frame(taxa = factor(final_order, levels = final_order))
constraint <- as.phylo(~taxa, data = constraint_data)
tree_reordered <- rotateConstr(tree, constraint)
tree_reordered <- ladderize(tree_reordered, right = FALSE)

# 绘制进化树
tree_plot <- ggtree(tree_reordered) + 
  theme_tree() +
  labs(title = paste("Tree (", length(final_order), " strains)"))

if (show_labels) {
  label_size <- max(3, 8 - 0.03 * length(final_order))
  tree_plot <- tree_plot + 
    geom_tiplab(size = label_size, align = TRUE, linesize = 0.2)
}
# 绘制基因存在缺失矩阵
# 准备矩阵绘图数据
plot_data <- data.frame(
  Gene = rownames(filtered_matrix),
  filtered_matrix,
  check.names = FALSE
) %>% melt(id.vars = "Gene", variable.name = "Strain", value.name = "Presence")
# 设置菌株顺序
plot_data$Strain <- factor(plot_data$Strain, levels = final_order)
plot_data$Gene <- factor(plot_data$Gene, levels = rownames(filtered_matrix))
# 根据菌株分组定义存在基因的颜色和对应的菌种名称
plot_data <- plot_data %>%
  mutate(
    Species = case_when(
      Presence == 0 ~ NA_character_,
      match(Strain, strain_order_268) <= index1 ~ "F.prausnitzii",
      match(Strain, strain_order_268) >= index2_start & 
        match(Strain, strain_order_268) <= index2_end ~ "F.duncaniae",
      match(Strain, strain_order_268) >= index3 ~ "F.longum"
    )
  )
# 绘制矩阵图
matrix_plot <- ggplot(plot_data, aes(x = Gene, y = Strain, fill = Species)) +
  geom_raster() +
  scale_fill_manual(
    values = setNames(species_colors$Color, species_colors$Species),
    na.value = "white",
    name = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = if(show_labels) element_text(size = label_size) else element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 10)
  ) +
  labs(title = paste("Genes Present in ≥", min_strains, "Strains\n(", 
                     nrow(filtered_matrix), " gene clusters)"))
# 组合图形
combined_plot <- plot_grid(
  tree_plot,
  matrix_plot,
  nrow = 1,
  align = "h",
  rel_widths = c(1, 3)
)
print(combined_plot)
