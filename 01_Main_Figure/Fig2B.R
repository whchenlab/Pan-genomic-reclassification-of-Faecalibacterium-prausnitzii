# Load required libraries and data
library(ggplot2)
library(minpack.lm)
matrix <- read.csv("gene_presence_absence.csv")

# Extract strain data (starting from column 4)
strains_data <- matrix[, 4:ncol(a)]
strain_names <- colnames(strains_data)
n_strains <- length(strain_names)
pan_core_data <- data.frame(
  n_strains = integer(),
  pan_genes = integer(),
  core_genes = integer()
)

# Set core gene threshold (90%)
CORE_THRESHOLD <- 0.90

# Calculate pan-genome and core-genome sizes for different numbers of strains
cat("Calculating pan-genome curve...\n")
for (n in 1:n_strains) {
  pan_values <- numeric(30)  # 30 replicate samples
  core_values <- numeric(30)
  
  for (i in 1:30) {
    selected_strains <- sample(strain_names, n)
    subset_data <- strains_data[, selected_strains, drop = FALSE]
    presence_matrix <- !apply(subset_data, 2, function(x) x == "")
    
    pan_values[i] <- sum(rowSums(presence_matrix) > 0)  # pan-genome size
    core_values[i] <- sum(rowSums(presence_matrix) >= (n * CORE_THRESHOLD))  # core-genome size
  }
  
  # Take average
  pan_core_data <- rbind(pan_core_data, data.frame(
    n_strains = n,
    pan_genes = mean(pan_values),
    core_genes = mean(core_values)
  ))
  
  # Show progress
  if (n %% 10 == 0 || n == n_strains) {
    cat(sprintf("Completed %d/%d (%.1f%%)\n", n, n_strains, n/n_strains*100))
  }
}

# Model fitting
cat("Fitting curves...\n")
# Pan-genome model: y = a*x^b + c
pan_model <- nlsLM(pan_genes ~ a * n_strains^b + c, 
                   data = pan_core_data,
                   start = list(a = 6000, b = 0.36, c = 3000),
                   control = nls.lm.control(maxiter = 500))

# Core-genome model: y = a*exp(b*x) + c
core_model <- nlsLM(core_genes ~ a * exp(b * n_strains) + c, 
                    data = pan_core_data,
                    start = list(a = 2500, b = -1, c = 1500),
                    control = nls.lm.control(maxiter = 500))

# Extract parameters
pan_params <- coef(pan_model)
core_params <- coef(core_model)

# Calculate R-squared
pan_pred <- predict(pan_model)
ss_res_pan <- sum((pan_core_data$pan_genes - pan_pred)^2)
ss_tot_pan <- sum((pan_core_data$pan_genes - mean(pan_core_data$pan_genes))^2)
r2_pan <- 1 - (ss_res_pan / ss_tot_pan)

core_pred <- predict(core_model)
ss_res_core <- sum((pan_core_data$core_genes - core_pred)^2)
ss_tot_core <- sum((pan_core_data$core_genes - mean(pan_core_data$core_genes))^2)
r2_core <- 1 - (ss_res_core / ss_tot_core)

# Create equation expressions
pan_eq <- paste0("Pan = ", round(pan_params["a"], 2), " * x^", round(pan_params["b"], 2), 
                 " + ", round(pan_params["c"], 2), ", R² = ", round(r2_pan, 3))

core_eq <- paste0("Core = ", round(core_params["a"], 2), " * exp(", round(core_params["b"], 3), 
                  " * x) + ", round(core_params["c"], 2), ", R² = ", round(r2_core, 3))

# Generate PDF and create plot
pdf_file <- "./pancurve.pdf"
p <- ggplot(pan_core_data, aes(x = n_strains)) +
  # Original data curves
  geom_line(aes(y = pan_genes, color = "Pan-genome"), linewidth = 1.2) +
  geom_line(aes(y = core_genes, color = "Core-genome"), linewidth = 1.2) +
  # Fitted curves
  stat_function(fun = function(x) pan_params["a"] * x^pan_params["b"] + pan_params["c"],
                aes(color = "Pan-genome fit"), linetype = 2, linewidth = 1) +
  stat_function(fun = function(x) core_params["a"] * exp(core_params["b"] * x) + core_params["c"],
                aes(color = "Core-genome fit"), linetype = 2, linewidth = 1) +
  # Equations
  annotate("text", x = 0.7 * max(pan_core_data$n_strains), 
           y = 0.9 * max(pan_core_data$pan_genes), 
           label = pan_eq, color = "blue", size = 5) +
  annotate("text", x = 0.7 * max(pan_core_data$n_strains), 
           y = 0.1 * max(pan_core_data$pan_genes), 
           label = core_eq, color = "red", size = 5) +
  # Core gene threshold annotation
  annotate("text", x = max(pan_core_data$n_strains) * 0.1, 
           y = max(pan_core_data$core_genes) * 0.5, 
           label = paste0("Core gene threshold: ", CORE_THRESHOLD*100, "%"),
           color = "darkgreen", size = 4, fontface = "bold") +
  
  labs(
    x = "Number of strains",
    y = "Number of genes",
    title = " ",
    subtitle = paste0("Core genes defined at ", CORE_THRESHOLD*100, "% presence threshold"),
    color = ""
  ) +
  scale_color_manual(values = c("Pan-genome" = "black", "Core-genome" = "black",
                                "Pan-genome fit" = "blue", "Core-genome fit" = "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# Save PDF
ggsave(pdf_file, p, width = 10, height = 8)
cat(sprintf("PDF saved: %s\n", pdf_file))
