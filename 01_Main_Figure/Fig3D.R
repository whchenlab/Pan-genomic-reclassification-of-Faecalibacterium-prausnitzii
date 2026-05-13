# Load required library and data
library(tidyverse)
data <- data.frame(
  Threshold = c("com70", "com75", "com80", "com85", "com90", "com95",
                "total_70", "total_75", "total_80", "total_85", "total_90", "total_95"),
  F.longum = c(1940, 1868, 1745, 1529, 1152, 588,
               2614, 2480, 2272, 1951, 1450, 680),
  F.prausnitzii = c(27, 27, 27, 26, 26, 18,
                    1776, 1678, 1559, 1410, 1149, 556),
  F.duncaniae = c(24, 23, 23, 23, 21, 16,
                  1448, 1343, 1208, 1085, 903, 479),
  F.taiwanense = c(33, 33, 33, 31, 27, 23,
                   1332, 1226, 1119, 1004, 858, 523),
  F.centralis = c(15, 14, 14, 13, 12, 8,
                  502, 420, 360, 271, 190, 59),
  F.indigenum = c(2, 2, 2, 2, 2, 2,
                  391, 363, 339, 311, 267, 140),
  F.hattorii = c(4, 4, 4, 4, 4, 3,
                 222, 217, 208, 198, 166, 83),
  F.bellum = c(1, 1, 1, 1, 1, 1,
               208, 159, 116, 90, 68, 35),
  F.butyricigenerans = c(110, 104, 101, 94, 70, 41,
                         134, 124, 120, 113, 91, 49),
  F.amidilyticum = c(0, 0, 0, 0, 0, 0,
                     67, 63, 61, 55, 49, 30),
  F.wellingii = c(1, 1, 1, 1, 1, 1,
                  67, 61, 58, 53, 44, 27),
  F.langellae = c(1, 1, 1, 1, 0, 0,
                  64, 62, 58, 55, 44, 16)
)

# Count Yes/No percent for each species
stack_data <- data.frame()
label_data <- data.frame()
strains <- colnames(data)[-1]

for(strain in strains) {
  # Extract data of >90 completeness genomes
  com90_value <- data[data$Threshold == "com90", strain]
  total_90_value <- data[data$Threshold == "total_90", strain]
  # Count percent
  Yes <- com90_value / total_90_value
  No <- 1 - Yes
  stack_data <- rbind(stack_data, 
                      data.frame(Strain = strain, Category = "Yes", Value = Yes),
                      data.frame(Strain = strain, Category = "No", Value = No))
  label_data <- rbind(label_data, 
                      data.frame(Strain = strain, Total = total_90_value))
}

# Show percent(>5%) on the bar plot
stack_data$Percentage_Label <- ifelse(
  stack_data$Value > 0.05,
  paste0(round(stack_data$Value * 100, 0), "%"),
  ""
)

# Set strain order
strain_order <- label_data %>% 
  arrange(desc(Total)) %>% 
  pull(Strain)

stack_data$Strain <- factor(stack_data$Strain, levels = strain_order)
label_data$Strain <- factor(label_data$Strain, levels = strain_order)

# Create bar plot
ggplot() +
  geom_bar(data = stack_data, 
           aes(x = Strain, y = Value, fill = Category),
           stat = "identity", width = 0.7) +
  
  # Add percent label on the bar
  geom_text(data = stack_data, 
            aes(x = Strain, y = Value, fill = Category, label = Percentage_Label),
            position = position_stack(vjust = 0.5),
            size = 2,
            color = "white",
            fontface = "bold") +
  
  # Add total strain number on the top
  geom_text(data = label_data, 
            aes(x = Strain, 
                y = 1.05,
                label = Total),
            vjust = 0,
            size = 3,
            color = "black",
            inherit.aes = FALSE) +

  scale_fill_manual(values = c("#FB9A99" ,"#E31A1C")) +

  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.15),
                     breaks = seq(0, 1, 0.2)) +
  
  labs(x = "Species", 
       y = "Proportion", 
       fill = "GalNac") +
  
  theme_minimal() +
  theme(
    text = element_text(family = "", size = 12, face = 'plain'),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0, 
                               family = '', size = 10, face = 'italic'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
ggsave(filename = "GalNac_com90_barplot.pdf", width = 6, height = 5)

