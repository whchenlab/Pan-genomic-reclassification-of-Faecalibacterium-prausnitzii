library(ggplot2)
library(dplyr)
library(scales)

# Merge Checkm and Busco datasets
combined_checkm_clean <- summary_checkm[summary_checkm$Completeness != "-", ]
combined_checkm_clean <- combined_checkm_clean %>%
  mutate(
    Completeness = as.numeric(Completeness),
    Contamination = as.numeric(Contamination),
    Method = "CheckM"
  )
combined_busco_clean <- summary_busco[summary_busco$Complete != "-", ]
combined_busco_clean <- combined_busco_clean %>%
  mutate(
    Complete = as.numeric(Complete),
    Fragmented = as.numeric(Fragmented),
    Duplicated = as.numeric(Duplicated),
    Method = "BUSCO"
  )

# Count filtered number by methods and thresholds
calculate_retained_genomes <- function(data, method_name) {
  total_genomes <- nrow(data)
  complete_thresholds <- seq(40, 100, by = 1)
  
  if (method_name == "CheckM") {
    # CheckM base standaed: Completeness >= threshold, Contamination <= 5
    retained_genomes <- sapply(complete_thresholds, function(threshold) {
      data %>%
        filter(Completeness >= threshold, Contamination <= 5) %>%
        nrow()
    })
  } else if (method_name == "BUSCO") {
    # BUSCO base standard: Complete >= threshold, Fragmented <= 10, Duplicated <= 5
    retained_genomes <- sapply(complete_thresholds, function(threshold) {
      data %>%
        filter(Complete >= threshold, Fragmented <= 10, Duplicated <= 5) %>%
        nrow()
    })
  }
  
  data.frame(
    Complete_Threshold = complete_thresholds,
    Retained_Genomes = retained_genomes,
    Method = method_name,
    Total_Genomes = total_genomes
  )
}

# Prepare plot data
checkm_data <- calculate_retained_genomes(combined_checkm_clean, "CheckM")
busco_data <- calculate_retained_genomes(combined_busco_clean, "BUSCO")
plot_data <- rbind(checkm_data, busco_data)

# Create plot
colors <- c("CheckM" = "#386CB0", "BUSCO" = "#7FC97F")
y_max <- max(plot_data$Retained_Genomes)
p <- ggplot(plot_data, aes(x = Complete_Threshold, y = Retained_Genomes, color = Method, group = Method)) +
  geom_line(linewidth = 1) +
  geom_point(data = plot_data[plot_data$Complete_Threshold %% 5 == 0, ], 
             size = 2) +
  geom_hline(data = plot_data %>% distinct(Method, Total_Genomes), 
             aes(yintercept = Total_Genomes, color = Method), 
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  geom_vline(xintercept = 50, color = "red", linetype = "solid", linewidth = 0.8) +
  geom_vline(xintercept = 70, color = "red", linetype = "solid", linewidth = 0.8) +
  geom_text(aes(x = 50, y = y_max * 0.3, label = "BUSCO Threshold (50)"), 
            color = "red", size = 3, fontface = "bold", vjust = 0, hjust = 0.5) +
  geom_text(aes(x = 70, y = y_max * 0.4, label = "CheckM Threshold (70)"), 
            color = "red", size = 3, fontface = "bold", vjust = 0, hjust = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "Faecalibacterium Genome Number with Different Complete Thresholds",
       x = "Complete Thresholds",
       y = "Remained Genomes Number",
       color = "Assessment Method",
       subtitle = "CheckM: Contamination ≤ 5 | BUSCO: Fragmented ≤ 10, Duplicated ≤ 5") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(40, 100))

for (method in unique(plot_data$Method)) {
  method_data <- plot_data[plot_data$Method == method & plot_data$Complete_Threshold %% 5 == 0, ]
  for (i in 1:nrow(method_data)) {
    p <- p + 
      geom_text(data = method_data[i, ], 
                aes(label = Retained_Genomes),
                vjust = -1, color = colors[method], size = 3, show.legend = FALSE)
  }
}

# Add supplement info
method_totals <- plot_data %>% 
  distinct(Method, Total_Genomes) %>%
  mutate(
    x_pos = 97,
    label = paste("Total:", Total_Genomes)
  )
p <- p + 
  geom_text(data = method_totals, 
            aes(x = x_pos, y = Total_Genomes, label = label, color = "black"),
            hjust = 0, size = 3.5, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = method_totals, 
            aes(x = 95.5, y = Total_Genomes, label = Method, color = Method),
            hjust = 1, size = 3.5, fontface = "bold",
            show.legend = FALSE)

# Print and save
print(p)
ggsave("Faecali_filter_number.pdf", 
       plot = p, 
       width = 8, 
       height = 6)
