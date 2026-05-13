# Load required package and data
library(tidyverse)
data <- data.frame(
  Threshold = c("com90", "total_90"),
  F.longum = c(1303, 1450),
  F.prausnitzii = c(1, 1149),
  F.duncaniae = c(1, 903),
  F.taiwanense = c(0, 858),
  F.centralis = c(0, 190),
  F.indigenum = c(0, 267),
  F.hattorii = c(0, 166),
  F.bellum = c(0, 68),
  F.butyricigenerans = c(87, 91),
  F.amidilyticum = c(0, 49),
  F.wellingii = c(0, 44),
  F.langellae = c(0, 43)
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
  if(total_90_value > 0) {
    Yes <- com90_value / total_90_value
  } else {
    Yes <- 0
  }
  No <- 1 - Yes
  stack_data <- rbind(stack_data, 
                      data.frame(Strain = strain, Category = "Yes", Value = Yes),
                      data.frame(Strain = strain, Category = "No", Value = No))
  label_data <- rbind(label_data, 
                      data.frame(Strain = strain, Total = total_90_value))
}

# Show percent(>5%) on the bar plot
stack_data$Percentage_Label <- ifelse(
  stack_data$Value > 0.05 & !is.na(stack_data$Value),
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
  
  scale_fill_manual(values = c("#FDBF6F", "#FF7F00")) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.15),
                     breaks = seq(0, 1, 0.2)) +
  
  labs(x = "Species", 
       y = "Proportion", 
       fill = "FAAH") +
  
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
ggsave(filename = "FAAH_com90_barplot.pdf", width = 6, height = 5)
