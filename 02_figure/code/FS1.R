library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

# 1. 导入excel
excel_path <- "meaphlan_SGB2GTDB_version.xlsx"  
sheet_names <- excel_sheets(excel_path)
excel_data <- lapply(sheet_names, function(sheet_name) {
  df <- read_excel(
    path = excel_path,
    sheet = sheet_name,
    col_names = FALSE
  )
  if (ncol(df) >= 1) colnames(df)[1] <- "SGB"       # 第一列命名为SGB
  if (ncol(df) >= 2) colnames(df)[2] <- "Taxonomy"  # 第二列命名为Taxonomy
  
  return(df)
})
names(excel_data) <- sheet_names

# 2. 预处理数据
processed_data <- lapply(excel_data, function(df) {
  df %>%
    mutate(SGB = str_remove(SGB, "_.*")) %>%
    mutate(Taxonomy = str_extract(Taxonomy, "(?<=;)[^;]+$"))
})
names(processed_data) <- names(excel_data)

# 3. 生成绘图表格
sheet_names <- names(processed_data)
new_colnames <- str_extract(sheet_names, "[^_]+$")
all_sgbs <- unlist(lapply(processed_data, function(df) df$SGB))
unique_sgbs <- unique(all_sgbs)
merged_data <- tibble(SGB = unique_sgbs)
# 按SGB匹配，依次合并5个表格的Taxonomy数据到merged_data
for (i in seq_along(processed_data)) {
  current_table <- processed_data[[i]] %>%
    distinct(SGB, .keep_all = TRUE) %>%
    select(SGB, Taxonomy) %>%
    rename(!!new_colnames[i] := Taxonomy)
  
  merged_data <- merged_data %>%
    left_join(current_table, by = "SGB")
}
# 新增ANI分类列
# (1) 创建SGB与ANI_Taxonomy的映射关系数据框
ani_mapping <- tibble(
  SGB = c(
    "SGB15315", "SGB15316", "SGB15318", "SGB15322", 
    "SGB15332", "SGB15340", "SGB15342", "SGB66208"
  ),
  ANI_Taxonomy = c(
    "F.butyricigenerans", "F.longum", "F.duncaniae", "F.hattorii",
    "F.prausnitzii", "F.wellingii", "F.faecis_F.taiwanense", "F.gallinarum"
  )
)
# (2) 合并映射关系到merged_data，并填充未匹配的行为"Unknown"
merged_data <- merged_data %>%
  left_join(ani_mapping, by = "SGB") %>%
  mutate(ANI_Taxonomy = replace_na(ANI_Taxonomy, "Unknown"))
# 数据清洗
merged_data <- merged_data %>%
  mutate(across(2:6, ~ str_replace(., "s__Faecalibacterium ", "F."))) %>%
  mutate(across(2:6, ~ str_replace(., "s__", "Undefined")))





library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)

# 一、数据处理与颜色映射准备
# 1. 定义颜色映射
prefix_colors <- c(
  "F.duncaniae" = "#F6FF33",        # 橙色
  "F.hattorii" = "#3357FF",         # 深蓝色
  "F.faecis_F.taiwanense" = "#33FFF6", # 紫色
  "F.longum" = "#FF3333",           # 红色
  "F.butyricigenerans" = "#33FF33", # 绿色
  "F.prausnitzii" = "#F633FF",      # 粉色
  "F.wellingii" = "#FFA500",        # 橙色
  "F.gallinarum" = "#8A2BE2",       # 靛蓝色
  "Unknown" = "#CCCCCC"             # 灰色
)
# 2. 为merged_data添加颜色列
merged_data <- merged_data %>%
  mutate(
    SGB_color = case_when(
      ANI_Taxonomy %in% names(prefix_colors) ~ prefix_colors[ANI_Taxonomy],
      TRUE ~ prefix_colors["Unknown"]
    )
  )

# 二、创建基础表格
# 1. 准备绘图数据（前6列）并将NA值替换为空字符串
plot_data <- merged_data %>% 
  select(1:6) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))
n_rows <- nrow(plot_data)  
# 2. 创建自定义主题
fill_colors <- matrix(merged_data$SGB_color, nrow = n_rows, ncol = ncol(plot_data))
original_data <- merged_data %>% select(1:6)
for(i in 1:nrow(original_data)) {
  for(j in 1:ncol(original_data)) {
    if(is.na(original_data[i, j])) {
      fill_colors[i, j] <- "white"
    }
  }
}

custom_theme <- ttheme_default(
  core = list(
    fg_params = list(fontsize = 9, hjust = 0.5),
    bg_params = list(
      fill = fill_colors,
      col = "black"
    )
  ),
  colhead = list(
    fg_params = list(fontsize = 10, fontface = "bold", hjust = 0.5),
    bg_params = list(fill = "#f0f0f0", col = "black")
  )
)
# 3. 生成基础table_grob
table_grob <- tableGrob(
  d = plot_data,
  theme = custom_theme,
  rows = NULL
)

# 三、添加"Metaphlan Version"标注
col_widths <- table_grob$widths
target_cols_width <- sum(col_widths[2:6])
target_cols_start <- col_widths[1]

version_label <- textGrob(
  label = "Metaphlan Version",
  x = unit(0.55, "npc"), 
  y = unit(0.95, "npc"),
  gp = gpar(fontsize = 15, fontface = "bold", col = "black"),
  just = "center"
)

table_with_label <- grobTree(
  table_grob,
  version_label,
  vp = viewport(width = 0.9, height = 0.95)
)

# 四、创建独立的图例
legend_data <- data.frame(
  Taxonomy = names(prefix_colors),
  Color = prefix_colors,
  stringsAsFactors = FALSE
)

create_legend <- function(legend_data) {
  n_items <- nrow(legend_data)
  legend_grob <- gTree(children = gList(
    rectGrob(
      x = unit(0, "npc"),
      y = unit(seq(0.8, 0.3, length.out = n_items), "npc"),
      width = unit(0.3, "npc"),
      height = unit(0.03, "npc"),
      gp = gpar(fill = legend_data$Color, col = "black")
    ),
    textGrob(
      label = legend_data$Taxonomy,
      x = unit(0.2, "npc"),
      y = unit(seq(0.8, 0.3, length.out = n_items), "npc"),
      just = "left",
      gp = gpar(fontsize = 9)
    ),
    textGrob(
      label = "ANI Cluster",
      x = unit(0.05, "npc"),
      y = unit(0.85, "npc"),
      just = "left",
      gp = gpar(fontsize = 12, fontface = "bold")
    )
  ))
  return(legend_grob)
}

legend_grob <- create_legend(legend_data)

# 五、组合表格和图例
final_with_legend <- arrangeGrob(
  table_with_label,
  legend_grob,
  nrow = 1,
  widths = c(0.85, 0.15),
  padding = unit(0, "lines")
)

# 六、在绘图窗口中显示结果
grid.newpage()
pushViewport(viewport(
  x = unit(0.5, "npc"),
  y = unit(0.5, "npc"),
  width = unit(1, "npc"),
  height = unit(1, "npc"),
  just = "center"
))
grid.draw(final_with_legend)
popViewport()
