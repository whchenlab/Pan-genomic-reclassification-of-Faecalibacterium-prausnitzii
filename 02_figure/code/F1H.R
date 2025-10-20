library(tidyverse)
library(ggplot2)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggdendro)

# --------------------------
# 1. 数据准备
# --------------------------
cog_description <- data.frame(
  Cat = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S", "X", "NA"),
  Description = c(
    "Translation, Ribosomal structure and biogenesis",
    "RNA processing and modification",
    "Transcription",
    "Replication, Recombination and Repair",
    "Chromatin structure and dynamics",
    "Cell cycle control. Cell division, Chromosome partitioning",
    "Nuclear structure",
    "Defense mechanisms",
    "Signal transduction mechanisms",
    "Cell wall/Membrane/Envelope biogenesis",
    "Cell motility",
    "Cytoskeleton",
    "Extracellular structures",
    "Intracellular Trafficking, Secretion,and Vesicular transport",
    "Posttranslational modification, Protein turnover, Chaperones",
    "Energy productionand conversion",
    "Carbohydrate transport and metabolism",
    "Amino acid transportand metabolism",
    "Nucleotide transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, Transfort and catabolism",
    "General function prepiction only",
    "Function unknown",
    "Nan",
    "Function unknown"
  ),
  BroadCategory = c(
    "Information Storage and Processing",  # J
    "Information Storage and Processing",  # A
    "Information Storage and Processing",  # K
    "Information Storage and Processing",  # L
    "Information Storage and Processing",  # B
    "Cellular Processes and Signaling",    # D
    "Cellular Processes and Signaling",    # Y
    "Cellular Processes and Signaling",    # V
    "Cellular Processes and Signaling",    # T
    "Cellular Processes and Signaling",    # M
    "Cellular Processes and Signaling",    # N
    "Cellular Processes and Signaling",    # Z
    "Cellular Processes and Signaling",    # W
    "Cellular Processes and Signaling",    # U
    "Cellular Processes and Signaling",    # O
    "Metabolism",                          # C
    "Metabolism",                          # G
    "Metabolism",                          # E
    "Metabolism",                          # F
    "Metabolism",                          # H
    "Metabolism",                          # I
    "Metabolism",                          # P
    "Metabolism",                          # Q
    "Information Storage and Processing",  # R
    "Function unknown",                    # S
    "Cellular Processes and Signaling",    # X
    "Function unknown"                     # NA
  )
)

broad_category_colors <- c(
  "Information Storage and Processing" = "#377EB8",
  "Cellular Processes and Signaling" = "#4DAF4A",
  "Metabolism" = "#E41A1C",
  "Function unknown" = "#999999"
)

# --------------------------
# 2. 核心组处理
# --------------------------
# 处理数据
eggnog_FP <- eggnog_FP %>% distinct(og_id, .keep_all = TRUE)  # 直接按og_id去重
eggnog_FD <- eggnog_FD %>% distinct(og_id, .keep_all = TRUE)
eggnog_FL <- eggnog_FL %>% distinct(og_id, .keep_all = TRUE)
eggnog_SHELL <- eggnog_SHELL %>% distinct(og_id, .keep_all = TRUE)
eggnog_CLOUD <- eggnog_CLOUD %>% distinct(og_id, .keep_all = TRUE)

# 定义基因集类型
common_core_genes <- eggnog_FP %>% filter(og_id %in% eggnog_FD$og_id & og_id %in% eggnog_FL$og_id) %>% mutate(Core_Type = "Common_Core")
common_shell_genes <- eggnog_SHELL %>% mutate(Core_Type = "Common_Shell")
common_cloud_genes <- eggnog_CLOUD %>% mutate(Core_Type = "Common_Cloud")
fp_core_genes <- eggnog_FP %>% mutate(Core_Type = "FP_Core")
fd_core_genes <- eggnog_FD %>% mutate(Core_Type = "FD_Core")
fl_core_genes <- eggnog_FL %>% mutate(Core_Type = "FL_Core")

# 合并所有基因集
all_core_genes <- bind_rows(
  common_core_genes, 
  common_shell_genes, 
  common_cloud_genes, 
  fp_core_genes, 
  fd_core_genes, 
  fl_core_genes
)

# --------------------------
# 3. 统计COG占比
# --------------------------
count_cog_categories <- function(df) {
  df %>%
    mutate(
      COG_category = ifelse(is.na(COG_category) | COG_category == "-", "NA", COG_category),
      COG_Letters = str_split(COG_category, "")
    ) %>%
    unnest(COG_Letters) %>%
    mutate(
      COG_Letters = ifelse(COG_Letters %in% tolower(cog_description$Cat), toupper(COG_Letters), COG_Letters),
      COG_Letters = ifelse(COG_Letters %in% cog_description$Cat, COG_Letters, "NA")
    ) %>%
    group_by(Core_Type, COG_Letters) %>%
    count() %>%
    ungroup() %>%
    complete(Core_Type, COG_Letters, fill = list(n = 0)) %>%
    left_join(cog_description, by = c("COG_Letters" = "Cat")) %>%
    rename(COG_Category = COG_Letters, Count = n)
}

# 重新计算COG计数
cog_counts <- count_cog_categories(all_core_genes)

# 计算每个基因集的总基因数
total_genes_per_core <- all_core_genes %>% 
  group_by(Core_Type) %>% 
  distinct(og_id) %>% 
  summarise(Total = n(), .groups = "drop")

# 计算比例，确保无缺失值
cog_proportions <- cog_counts %>%
  left_join(total_genes_per_core, by = "Core_Type") %>%
  mutate(
    Proportion = Count / Total * 100,
    Proportion = ifelse(is.na(Proportion), 0, Proportion)
  ) %>%
  select(Core_Type, Description, BroadCategory, Proportion)

# 设置x轴顺序
core_order <- c("Common_Core", "Common_Shell", "Common_Cloud", "FP_Core", "FD_Core", "FL_Core")
cog_proportions$Core_Type <- factor(cog_proportions$Core_Type, levels = core_order)

broad_category_summary <- cog_proportions %>%
  group_by(Core_Type, BroadCategory) %>%
  summarise(TotalProportion = sum(Proportion), .groups = "drop")

# --------------------------
# 4. 绘图元素准备
# --------------------------
# 转换为宽格式用于聚类
heatmap_data_wide <- cog_proportions %>%
  select(Core_Type, Description, Proportion) %>%
  pivot_wider(names_from = Core_Type, values_from = Proportion, values_fill = 0) %>%
  as.data.frame()

rownames(heatmap_data_wide) <- heatmap_data_wide$Description
heatmap_data_wide$Description <- NULL
heatmap_matrix <- as.matrix(heatmap_data_wide)

# 行聚类（功能类别）
row_clust <- hclust(dist(heatmap_matrix, method = "euclidean"), method = "average")
row_order <- rownames(heatmap_matrix)[row_clust$order]

# 调整因子水平
cog_proportions$Description <- factor(cog_proportions$Description, levels = row_order)
color_data <- cog_proportions %>% distinct(Description, BroadCategory) 
color_data$Description <- factor(color_data$Description, levels = row_order)

# 创建聚类树grob
dend_data <- ggdendro::dendro_data(row_clust)
dend_plot <- ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = unit(c(0,0,0,0), "mm"))
dend_grob <- ggplotGrob(dend_plot)

# 热图（0.0%-50.0%颜色图例）
heatmap_plot <- ggplot(cog_proportions, aes(x = Core_Type, y = Description, fill = Proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), size = 3, color = "black") +
  scale_fill_gradient(
    low = "#E6F2FF", 
    high = "#1E3A8A", 
    na.value = "white", 
    name = "Percentage",
    limits = c(0, 50),
    breaks = seq(0, 50, by = 10)
  ) +
  scale_y_discrete(limits = levels(cog_proportions$Description)) +
  labs(x = "Gene Set", y = "") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.margin = unit(c(0,0,0,0), "mm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 5)
  )

# 提取热图颜色图例
heatmap_legend <- get_legend(heatmap_plot)
heatmap_plot <- heatmap_plot + theme(legend.position = "none")

# 左侧颜色块
category_legend <- ggplot(color_data, aes(x = 1, y = Description, fill = BroadCategory)) +
  geom_tile(color = "white", linewidth = 0.5, width = 0.8) +
  scale_fill_manual(values = broad_category_colors) +
  scale_y_discrete(limits = levels(color_data$Description)) +
  theme_void() + 
  theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm"))

# 右侧标注（功能描述）
desc_annot_data <- color_data %>%
  mutate(y_pos = as.numeric(Description)) %>%
  arrange(y_pos)

right_annot_grob <- grobTree(
  textGrob(
    label = desc_annot_data$Description,
    x = 0, y = desc_annot_data$y_pos/max(desc_annot_data$y_pos),
    just = "left", gp = gpar(fontsize = 7.5),
    hjust = 0
  )
)

# 饼图
plot_pie <- function(core_type) {
  data <- broad_category_summary %>% filter(Core_Type == core_type)
  ggplotGrob(
    ggplot(data, aes(x = "", y = TotalProportion, fill = BroadCategory)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = broad_category_colors) +
      theme_void() + theme(legend.position = "none")
  )
}
pie_list <- map(core_order, plot_pie)

# 标题
main_title <- textGrob("COG Functional Distribution of Core/Shell/Cloud Genes",
                       gp = gpar(fontface = "bold", fontsize = 14))

# --------------------------
# 5. 调整位置参数
# --------------------------
dend_x <- 0.15        # 聚类树x坐标
dend_width <- 0.15    # 聚类树宽度
dend_ylow <- 0.182    # 聚类树底部y坐标
dend_yhigh <- 0.86    # 聚类树顶部y坐标

heatmap_x <- 0.45     # 热图x坐标
heatmap_y <- 0.5
heatmap_width <- 0.35 # 热图宽度
heatmap_height <- 0.7

color_block_x <- 0.25  # 颜色块x坐标
color_block_ylim <- c(0.185, 0.85)

right_annot_x <- 0.8  # 右侧标注x坐标
right_annot_width <- 0.3  # 右侧标注宽度

pie_x <- c(0.31, 0.37, 0.43, 0.49, 0.55, 0.61)  # 饼图x坐标
pie_y <- 0.91
pie_size <- 0.06

legend_x <- 0.06  
legend_start_y <- 0.95
legend_gap <- 0.03

# 热图颜色图例参数
heatmap_legend_x <- 0.1    # 图例x坐标
heatmap_legend_y <- 0.15    # 图例y坐标
heatmap_legend_width <- 0.08  # 图例宽度
heatmap_legend_height <- 0.6  # 图例高度

# --------------------------
# 6. 绘制所有元素
# --------------------------
grid.newpage()

# 绘制标题
pushViewport(viewport(x = 0.5, y = 0.1, height = 0.05))
grid.draw(main_title)
upViewport()

# 绘制6个饼图
for (i in 1:6) {
  pushViewport(viewport(
    x = pie_x[i], 
    y = pie_y, 
    width = pie_size, 
    height = pie_size
  ))
  grid.draw(pie_list[[i]])
  upViewport()
}

# 绘制聚类树
pushViewport(viewport(
  x = dend_x, 
  y = (dend_ylow + dend_yhigh)/2, 
  width = dend_width, 
  height = dend_yhigh - dend_ylow,
  just = "center"
))
grid.draw(dend_grob)
upViewport()

# 绘制左侧颜色块
pushViewport(viewport(
  x = color_block_x, 
  y = mean(color_block_ylim), 
  width = 0.05, 
  height = diff(color_block_ylim)
))
grid.draw(ggplotGrob(category_legend))
upViewport()

# 绘制热图
pushViewport(viewport(
  x = heatmap_x, 
  y = heatmap_y, 
  width = heatmap_width, 
  height = heatmap_height
))
grid.draw(ggplotGrob(heatmap_plot))
upViewport()

# 绘制右侧标注
pushViewport(viewport(
  x = right_annot_x, 
  y = 0.5, 
  width = right_annot_width, 
  height = 0.66,
  just = "center"
))
grid.draw(right_annot_grob)
upViewport()

# 绘制Broad Category图例
grid.text(
  "Broad Category",
  x = legend_x,
  y = legend_start_y + legend_gap,
  gp = gpar(fontface = "bold", fontsize = 10),
  just = "left"
)

for (i in seq_along(broad_category_colors)) {
  current_y <- legend_start_y - (i - 1) * legend_gap
  grid.rect(
    x = unit(legend_x, "npc"),
    y = unit(current_y, "npc"),
    width = unit(8, "mm"),
    height = unit(4, "mm"),
    gp = gpar(fill = broad_category_colors[i], col = "black", lwd = 1),
    just = "left"
  )
  grid.text(
    names(broad_category_colors)[i],
    x = unit(legend_x + 0.04, "npc"),
    y = unit(current_y, "npc"),
    gp = gpar(fontsize = 9),
    just = "left"
  )
}

# 绘制热图颜色图例（0.0%-50.0%）
pushViewport(viewport(
  x = heatmap_legend_x,
  y = heatmap_legend_y,
  width = heatmap_legend_width,
  height = heatmap_legend_height,
  just = "center"
))
grid.draw(heatmap_legend)
upViewport()

