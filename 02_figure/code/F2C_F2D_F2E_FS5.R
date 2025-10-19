library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(cowplot)

# 定义第一类CSV处理函数（基础文件：un19和un11）
generate_plot <- function(file_path, strain_title, is_first) {
  data <- read.csv(file_path, skip = 1, header = TRUE)
  
  # 菌株名处理
  strain_names <- data %>% 
    select(1) %>%
    pull() %>%
    make.unique(sep = "_")
  n_strains <- length(strain_names)
  
  # 目标基因列
  target_genes <- colnames(data)[-1]
  
  # 转换为数值型
  data_numeric <- data %>%
    select(-1) %>%
    mutate(across(everything(), as.numeric))
  
  # 颜色分类逻辑
  strain_colors <- character(nrow(data_numeric))
  for (i in 1:nrow(data_numeric)) {
    row_data <- data_numeric[i, ]
    cols_2535_2545 <- row_data[7:17]
    na_count <- sum(is.na(cols_2535_2545))
    
    if (na_count >= 9) {
      strain_colors[i] <- '#FD2598'
    } else if (is.na(row_data[1]) && is.na(row_data[2])) {
      strain_colors[i] <- '#1D73FF'
    } else {
      vals <- unlist(row_data[1:17])[!is.na(unlist(row_data[1:17]))]
      if (length(vals) > 1) {
        diffs <- diff(vals)
        strain_colors[i] <- ifelse(all(diffs >= 0) | all(diffs <= 0), '#7FC97F', '#FDD80A')
      } else {
        strain_colors[i] <- '#FDD80A'
      }
    }
  }
  
  # 绘图数据框
  plot_data <- data_numeric %>%
    mutate(Strain = strain_names, Color = strain_colors) %>%
    pivot_longer(cols = -c(Strain, Color), names_to = "Gene", values_to = "Value") %>%
    mutate(Present = !is.na(Value), Color = ifelse(Present, Color, "white"))
  
  # 颜色标签
  color_labels <- c(
    "white" = "Not Available",
    '#FD2598' = "GalNac Cluster Absent",
    '#1D73FF' = "Insertion Direction Unknown",
    '#7FC97F' = "GalNac Forward Insertion",
    '#FDD80A' = "GalNac Reverse Insertion"
  )
  plot_data$Color <- factor(plot_data$Color, levels = names(color_labels), labels = color_labels)
  
  # 菌株聚类
  cluster_data <- data_numeric[, 1:17]
  cluster_data[is.na(cluster_data)] <- -10000
  dist_matrix <- dist(cluster_data, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")
  ordered_strains <- strain_names[hc$order]
  plot_data$Strain <- factor(plot_data$Strain, levels = ordered_strains)
  
  # 基因标签
  desired_genes <- c(
    "dinB", "2530", "lacC", "2532", "2533", "ptsH", "gatZ-kbaZ", "nagA", 
    "agaS", "agaF", "2539", "gatY-kbaY", "GH109", "rhaR", "agaV", "agaC", 
    "agaD", "2546", "2547", "2548", "2549", "2550", "2551", "2552", 
    "2553", "2554", "2555", "2556", "2557", "2558", "2559", "2560"
  )
  gene_labels <- ifelse(grepl("^[a-zA-Z]", desired_genes), desired_genes, "")
  plot_data$Gene <- factor(plot_data$Gene, levels = target_genes)
  
  # 判断是否为基础文件组最后一个图
  is_last_base <- strain_title == "Unknown_11"
  
  # 绘图（删除基础组x轴标题）
  p <- ggplot(plot_data, aes(x = Gene, y = Strain, fill = Color)) +
    geom_tile(color = "gray30", linewidth = 0.2) +
    scale_fill_manual(values = c(
      "Not Available" = "white",
      "GalNac Cluster Absent" = '#FD2598',
      "Insertion Direction Unknown" = '#1D73FF',
      "GalNac Forward Insertion" = '#7FC97F',
      "GalNac Reverse Insertion" = '#FDD80A'
    ), drop = FALSE, name = "Legend") +
    scale_x_discrete(breaks = target_genes, labels = gene_labels) +
    labs(
      y = strain_title,
      x = if (is_last_base) "" else ""  # 删除基础组x轴标题
    ) +
    theme_minimal(base_size = 10) +
    theme(
      text = element_text(family = "SimHei"),
      axis.text.x = if (is_last_base) {
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)  # 显示横坐标标签
      } else {
        element_blank()
      },
      axis.ticks.x = if (is_last_base) element_line() else element_blank(),
      axis.title.x = if (is_last_base) element_blank() else element_blank(),  # 不显示x轴标题
      axis.text.y = element_text(size = 4, angle = 0, hjust = 1),
      axis.title.y = element_text(size = 7, angle = 0, vjust = 0.5, hjust = 1),
      legend.position = if (is_first) "bottom" else "none",  # 基础组第一个图显示图例
      panel.grid = element_blank(),
      plot.margin = margin(0, 5.5, 0, 5.5, "pt")
    )
  
  return(list(plot = p, n = n_strains))
}

# 定义第二类CSV处理函数（特殊文件：un6和un12）
generate_special_plot <- function(file_path, strain_title, is_first_special) {
  # 目标基因顺序
  desired_genes <- c(
    "1666", "RspR", "1668", "1669", "1670", "ptsH", "gatZ-kbaZ", "nagA", 
    "agaS", "agaF", "1676", "gatY-kbaY", "GH109", "rhaR", "agaV", "agaC", 
    "agaD", "1683", "1684", "XerC", "1686", "1687", "1688", "1689", "YwiE", 
    "1691", "1692", "1693", "1694"
  )
  
  # 基因标签
  gene_labels <- ifelse(grepl("^[a-zA-Z]", desired_genes), desired_genes, "")
  
  # 读取数据
  data <- read.csv(file_path, skip = 1, header = TRUE)
  
  # 菌株名提取
  strain_names <- data %>% 
    select(1) %>%
    pull() %>%
    make.unique(sep = "_")
  n_strains <- length(strain_names)
  
  # 分离菌株列和基因数据列
  gene_data <- data %>% select(-1)
  
  # 确保基因列名与desired_genes一致
  if (ncol(gene_data) == length(desired_genes)) {
    colnames(gene_data) <- desired_genes
  } else {
    warning("基因列数量与desired_genes不一致！请检查数据。")
    desired_genes <- colnames(gene_data)
    gene_labels <- ifelse(grepl("^[a-zA-Z]", desired_genes), desired_genes, "")
  }
  
  # 转换基因数据为数值型
  data_numeric <- gene_data %>%
    mutate(across(everything(), as.numeric))
  
  # 颜色分类逻辑
  strain_colors <- character(nrow(data_numeric))
  for (i in 1:nrow(data_numeric)) {
    row_data <- data_numeric[i, ]
    cols_2535_2545 <- row_data[7:17]
    na_count <- sum(is.na(cols_2535_2545))
    
    if (na_count >= 9) {
      strain_colors[i] <- '#FD2598'
    } else if (is.na(row_data[1]) && is.na(row_data[2])) {
      strain_colors[i] <- '#1D73FF'
    } else {
      vals <- unlist(row_data[1:17])[!is.na(unlist(row_data[1:17]))]
      if (length(vals) > 1) {
        diffs <- diff(vals)
        strain_colors[i] <- ifelse(all(diffs >= 0) | all(diffs <= 0), '#7FC97F', '#FDD80A')
      } else {
        strain_colors[i] <- '#FDD80A'
      }
    }
  }
  
  # 创建绘图数据框
  plot_data <- data_numeric %>%
    mutate(Color = strain_colors) %>%
    pivot_longer(
      cols = -Color,
      names_to = "Gene",
      values_to = "Value",
      values_drop_na = FALSE
    ) %>%
    mutate(
      Strain = rep(strain_names, each = length(desired_genes)),
      Gene = factor(Gene, levels = desired_genes),
      Present = !is.na(Value),
      Color = ifelse(Present, Color, "white")
    )
  
  # 颜色标签
  color_labels <- c(
    "white" = "Not Available",
    '#FD2598' = "GalNac Cluster Absent",
    '#1D73FF' = "Insertion Direction Unknown",
    '#7FC97F' = "GalNac Forward Insertion",
    '#FDD80A' = "GalNac Reverse Insertion"
  )
  plot_data$Color <- factor(plot_data$Color, levels = names(color_labels), labels = color_labels)
  
  # 菌株聚类
  cluster_data <- data_numeric[, 1:17]
  cluster_data[is.na(cluster_data)] <- -10000
  dist_matrix <- dist(cluster_data, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")
  ordered_strains <- strain_names[hc$order]
  plot_data$Strain <- factor(plot_data$Strain, levels = ordered_strains)
  
  # 判断是否为特殊文件组最后一个图
  is_last_special <- strain_title == "Unknown_12"
  
  # 绘图（修改特殊组x轴标题为"Target Genes"）
  p <- ggplot(plot_data, aes(x = Gene, y = Strain, fill = Color)) +
    geom_tile(color = "gray30", linewidth = 0.2) +
    scale_fill_manual(
      values = c(
        "Not Available" = "white",
        "GalNac Cluster Absent" = '#FD2598',
        "Insertion Direction Unknown" = '#1D73FF',
        "GalNac Forward Insertion" = '#7FC97F',
        "GalNac Reverse Insertion" = '#FDD80A'
      ),
      drop = FALSE,
      name = "Legend"
    ) +
    scale_x_discrete(
      breaks = desired_genes,
      labels = gene_labels
    ) +
    labs(
      y = strain_title,
      x = if (is_last_special) "Target Genes" else ""  # 修改为Target Genes
    ) +
    theme_minimal(base_size = 10) +
    theme(
      text = element_text(family = "SimHei"),
      axis.text.x = if (is_last_special) {
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
      } else {
        element_blank()
      },
      axis.text.y = element_text(size = 4, angle = 0, hjust = 1),
      axis.title.y = element_text(size = 7, angle = 0, vjust = 0.5, hjust = 1),
      axis.ticks.x = if (is_last_special) element_line() else element_blank(),
      axis.title.x = if (is_last_special) element_text(size = 10) else element_blank(),
      legend.position = if (is_first_special) "bottom" else "none",  # 特殊组第一个图显示图例
      panel.grid = element_blank(),
      plot.margin = margin(0, 5.5, 0, 5.5, "pt")
    )
  
  return(list(plot = p, n = n_strains))
}

# ----------------------
# 第一部分：处理基础CSV文件（un19和un11）
# ----------------------
files <- c("output_un19_1.csv", "output_un11_1.csv")
titles <- c("Unknown_19", "Unknown_11")
is_first_flags <- c(TRUE, FALSE)  # 基础组第一个图显示图例

base_results <- mapply(
  generate_plot, 
  file_path = files, 
  strain_title = titles, 
  is_first = is_first_flags, 
  SIMPLIFY = FALSE
)
base_plots <- lapply(base_results, function(x) x$plot)
n_base_strains <- sapply(base_results, function(x) x$n)

# 提取基础组图例
base_legend <- get_legend(base_plots[[1]] + guides(fill = guide_legend(nrow = 1)))
base_plots[[1]] <- base_plots[[1]] + theme(legend.position = "none")


# ----------------------
# 第二部分：处理特殊CSV文件（un6和un12）
# ----------------------
special_files <- c("output_un6.csv", "output_un12.csv")
special_titles <- c("Unknown_6", "Unknown_12")
is_first_special_flags <- c(TRUE, FALSE)  # 特殊组第一个图显示图例

special_results <- mapply(
  generate_special_plot, 
  file_path = special_files, 
  strain_title = special_titles, 
  is_first_special = is_first_special_flags, 
  SIMPLIFY = FALSE
)
special_plots <- lapply(special_results, function(x) x$plot)
n_special_strains <- sapply(special_results, function(x) x$n)

# 提取特殊组图例
special_legend <- get_legend(special_plots[[1]] + guides(fill = guide_legend(nrow = 1)))
special_plots[[1]] <- special_plots[[1]] + theme(legend.position = "none")


# ----------------------
# 组合所有图形
# ----------------------
all_plots <- c(base_plots, special_plots)

# 计算高度比例
total_strains <- sum(n_base_strains) + sum(n_special_strains)
plot_heights <- c(n_base_strains, n_special_strains) / total_strains * 0.8
legend_height <- 0.1
all_heights <- c(plot_heights, legend_height, legend_height)

# 组合图形
final_plot <- wrap_plots(
  c(all_plots, list(base_legend, special_legend)),
  ncol = 1, 
  heights = all_heights
) & theme(plot.margin = margin(0, 0, 0, 0, "pt"))

# 显示并保存图形
print(final_plot)
# ggsave("final_combined_plot.png", final_plot, width = 14, height = 12, dpi = 300)
