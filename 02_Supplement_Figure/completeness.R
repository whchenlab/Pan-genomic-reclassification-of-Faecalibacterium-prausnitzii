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

# Prepare plot data
checkm_plot_data <- combined_checkm_clean %>%
  select(Completeness, Method) %>%
  rename(Completeness_Value = Completeness)

busco_plot_data <- combined_busco_clean %>%
  select(Complete, Method) %>%
  rename(Completeness_Value = Complete)
plot_data <- rbind(checkm_plot_data, busco_plot_data)

# Create density plot
colors <- c("CheckM" = "#386CB0", "BUSCO" = "#7FC97F")
p_density <- ggplot(plot_data, aes(x = Completeness_Value, fill = Method)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = colors) +
  labs(title = "",
       x = "Completeness (%)",
       y = "Density",
       fill = "Assessment Method") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5))

# Add supplement info
method_stats <- plot_data %>%
  group_by(Method) %>%
  summarise(
    Total_Genomes = n(),
    Mean_Completeness = round(mean(Completeness_Value), 2),
    Median_Completeness = round(median(Completeness_Value), 2),
    label = paste0("Total: ", Total_Genomes, "\n",
                   "Mean: ", Mean_Completeness, "%\n",
                   "Median: ", Median_Completeness, "%")
  )
p_density <- p_density + 
  geom_text(data = method_stats, 
            aes(x = 85, y = c(Inf, Inf), 
                label = label, 
                color = Method),
            vjust = c(2, 4), hjust = 0, size = 2.5, fontface = "bold",
            show.legend = FALSE) +
  scale_color_manual(values = c( "#3F007D","#00441B"))

# Print and save
print(p_density)
ggsave("Faecali_completeness.pdf", 
       plot = p_density, 
       width = 8, 
       height = 6)
