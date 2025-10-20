library(ggplot2)
library(minpack.lm)
a <- read.csv("gene_presence_absence.csv")

# 提取菌株数据（从第15列）
strains_data <- a[, 15:ncol(a)]
strain_names <- colnames(strains_data)
n_strains <- length(strain_names)
pan_core_data <- data.frame(
  n_strains = integer(),
  pan_genes = integer(),
  core_genes = integer()
)

# 计算不同菌株数量下的泛基因和核心基因数量
for (n in 1:n_strains) {
  pan_values <- numeric(50)  # 50次重复抽样
  core_values <- numeric(50)
  
  for (i in 1:50) {
    selected_strains <- sample(strain_names, n)
    subset_data <- strains_data[, selected_strains, drop = FALSE]
    presence_matrix <- !apply(subset_data, 2, function(x) x == "")  # 缺失
    
    pan_values[i] <- sum(rowSums(presence_matrix) > 0)  # 泛基因数
    core_values[i] <- sum(rowSums(presence_matrix) == n)  # 核心基因数
  }
  
  # 取平均值
  pan_core_data <- rbind(pan_core_data, data.frame(
    n_strains = n,
    pan_genes = mean(pan_values),
    core_genes = mean(core_values)
  ))
}

# 计算拟合
# 泛基因拟合 - 幂函数模型: y = a*x^b + c
pan_model <- nlsLM(pan_genes ~ a * n_strains^b + c, 
                   data = pan_core_data,
                   start = list(a = 100, b = 0.5, c = 1000))
# 核心基因拟合 - 指数衰减模型: y = a*exp(b*x) + c
core_model <- nlsLM(core_genes ~ a * exp(b * n_strains) + c, 
                    data = pan_core_data,
                    start = list(a = -1000, b = -0.1, c = 3000))
# 提取拟合参数
pan_params <- coef(pan_model)
core_params <- coef(core_model)

# 计算R²
# 泛基因组模型R²
pan_pred <- predict(pan_model)
ss_res_pan <- sum((pan_core_data$pan_genes - pan_pred)^2)
ss_tot_pan <- sum((pan_core_data$pan_genes - mean(pan_core_data$pan_genes))^2)
r2_pan <- 1 - (ss_res_pan / ss_tot_pan)
# 核心基因组模型R²
core_pred <- predict(core_model)
ss_res_core <- sum((pan_core_data$core_genes - core_pred)^2)
ss_tot_core <- sum((pan_core_data$core_genes - mean(pan_core_data$core_genes))^2)
r2_core <- 1 - (ss_res_core / ss_tot_core)

# 构建表达式
pan_eq <- substitute(italic(Pan) == a * italic(x)^b + c * "," ~~ italic(R^2) == r2,
                     list(a = format(pan_params["a"], digits = 3),
                          b = format(pan_params["b"], digits = 3),
                          c = format(pan_params["c"], digits = 3),
                          r2 = format(r2_pan, digits = 3)))
pan_eq <- as.expression(pan_eq)

core_eq <- substitute(italic(Core) == a * e^(b * italic(x)) + c * "," ~~ italic(R^2) == r2,
                      list(a = format(core_params["a"], digits = 3),
                           b = format(core_params["b"], digits = 3),
                           c = format(core_params["c"], digits = 3),
                           r2 = format(r2_core, digits = 3)))
core_eq <- as.expression(core_eq)

# 绘制图形
ggplot(pan_core_data, aes(x = n_strains)) +
  # 原始数据曲线
  geom_line(aes(y = pan_genes, color = "Pan-genome"), linewidth = 1.2) +
  geom_line(aes(y = core_genes, color = "Core-genome"), linewidth = 1.2) +
  # 拟合曲线
  stat_function(fun = function(x) pan_params["a"] * x^pan_params["b"] + pan_params["c"],
                aes(color = "Pan-genome fit"), linetype = 2, linewidth = 1) +
  stat_function(fun = function(x) core_params["a"] * exp(core_params["b"] * x) + core_params["c"],
                aes(color = "Core-genome fit"), linetype = 2, linewidth = 1) +
  # 表达式
  annotate("text", x = 0.7 * max(pan_core_data$n_strains), 
           y = 0.9 * max(pan_core_data$pan_genes), 
           label = pan_eq, parse = TRUE, color = "blue", size = 5) +
  annotate("text", x = 0.7 * max(pan_core_data$n_strains), 
           y = 0.1 * max(pan_core_data$pan_genes), 
           label = core_eq, parse = TRUE, color = "red", size = 5) +
  
  labs(
    x = "Number of strains",
    y = "Number of genes",
    title = "Pangenome Curve Analysis of F.prausnitzii",
    color = ""
  ) +
  scale_color_manual(values = c("Pan-genome" = "black", "Core-genome" = "black",
                                "Pan-genome fit" = "blue", "Core-genome fit" = "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

