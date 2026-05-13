library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# ==================== 1. Import data ====================
# Species, Total, Thla, Hbd, Cro, Bcd, But
gene_data_raw <- data.frame(
  Species = c("F.butyricigenerans", "F.duncaniae", "F.faecis", "F.hattorii", 
              "F.langellae", "F.longum", "F.prausnitzii", "F.centralis", 
              "F.indigenum", "F.bellum", "F.amidilyticum", "F.wellingii"),
  Total   = c(91, 903, 858, 166, 44, 1450, 1149, 190, 267, 68, 49, 44),
  Thla    = c(75, 729, 696, 142, 29, 1217, 875, 147, 214, 54, 39, 33),
  Hbd     = c(65, 576, 725, 146, 29, 1203, 807, 133, 203, 52, 45, 37),
  Cro     = c(72, 704, 742, 154, 32, 1228, 881, 140, 237, 55, 44, 34),
  Bcd     = c(62, 561, 668, 130, 31, 1088, 726, 116, 202, 56, 40, 35),
  But     = c(85, 859, 840, 156, 41, 1392, 1046, 168, 246, 65, 45, 40)
)

# Define species name
all_strains <- c("F.butyricigenerans", "F.duncaniae", "F.faecis", "F.hattorii",
                 "F.langellae", "F.longum", "F.prausnitzii", "F.centralis",
                 "F.indigenum", "F.bellum", "F.amidilyticum", "F.wellingii")

# ==================== 2. Convert format ====================
# Define gene name
gene_long <- gene_data_raw %>%
  pivot_longer(
    cols = c(Thla, Hbd, Cro, Bcd, But),
    names_to = "Gene",
    values_to = "Yes"
  ) %>%
  mutate(
    Gene = recode(Gene, Thla = "thlA"),
    Gene = recode(Gene, Hbd = "hbd"),
    Gene = recode(Gene, Cro = "crt"),
    Gene = recode(Gene, Bcd = "bcd"),
    Gene = recode(Gene, But = "but"),
    Gene = factor(Gene, levels = c("thlA", "hbd", "crt", "bcd", "but"))
  )

# Extract total number
total_90 <- gene_data_raw %>%
  select(Species, Total) %>%
  rename(Group = Species, Total_Number = Total)

# Count Yes/No percent
sum_df_new <- gene_long %>%
  left_join(total_90, by = c("Species" = "Group")) %>%
  mutate(
    Yes = as.numeric(Yes),
    Total = as.numeric(Total_Number),
    No = Total - Yes,
    Yes_percent = Yes / Total * 100,
    No_percent = No / Total * 100
  ) %>%
  select(Gene, Species, Yes_percent, No_percent) %>%
  pivot_longer(
    cols = c(Yes_percent, No_percent),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    Category = ifelse(Category == "Yes_percent", "Yes", "No"),
    Species = factor(Species, levels = rev(all_strains))
  )

# ==================== 3. Generate heatmap ====================
custom_colors <- c("#CAB2D6", "#6A3D9A")

stacked_barplot_sorted <- ggplot(sum_df_new, aes(x = Percentage, y = Species, fill = Category)) +
  geom_bar(stat = "identity", width = 0.25) +
  facet_wrap(~ Gene, nrow = 1) +
  theme_minimal() +
  labs(x = "Percentage (%)", y = "Species", fill = "Category") +
  theme(
    text = element_text(size = 11),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 9, face = "italic"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_rect(fill = "grey95", color = NA)
  ) +
  geom_text(
    aes(label = ifelse(Percentage > 5, paste0(round(Percentage, 0), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))

print(stacked_barplot_sorted)

# Save as pdf
ggsave(filename = "Butyrate_com90_barplot_pangenome.pdf", 
       plot = stacked_barplot_sorted, 
       width = 8, 
       height = 6,
       units = "in")
