library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Import data
gene_data_raw <- read.table(text = "
Gene	F.longum	F.prausnitzii	F.duncaniae	F.faecis	F.un1	F.un2	F.hattorii	F.un3	F.butyricigenerans	F.un4	F.wellingii	F.langellae
bcd_90 threshold	1704	973	803	873	248	286	159	154	87	49	47	41
bcd_80 threshold	1708	979	808	874	248	286	159	154	87	50	47	41
bcd_70 threshold	1724	1030	820	882	250	305	170	161	87	52	47	43
bcd_60 threshold	1721	1021	819	878	250	305	171	159	87	51	47	43
bcd_50 threshold	1790	1080	843	886	251	311	179	159	95	51	47	43
bcd_COG	1765	1022	865	735	256	272	179	162	93	52	46	43
but_90 threshold	1044	1497	1316	1240	411	343	206	188	120	0	56	57
but_80 threshold	2396	1515	1327	1250	415	347	207	190	124	58	57	58
but_70 threshold	2441	1544	1348	1260	418	349	210	192	126	59	57	59
but_60 threshold	2408	1518	1334	1251	416	348	208	190	124	58	57	59
but_50 threshold	2408	1518	1334	1253	416	348	208	190	124	59	57	59
but_COG	2414	1593	1331	1276	418	355	215	192	125	58	55	56
cro_90 threshold	1974	1203	1024	989	308	323	198	157	102	55	49	40
cro_80 threshold	1982	1210	1032	994	309	330	198	158	102	55	49	40
cro_70 threshold	1993	1217	1039	1004	312	330	198	158	102	55	50	40
cro_60 threshold	1989	1214	1035	997	309	330	198	158	102	55	50	40
cro_50 threshold	1989	1214	1037	998	309	330	198	158	102	55	50	40
cro_COG	1934	1209	1040	878	263	323	179	153	104	54	50	39
hbd_90 threshold	1913	1089	834	954	277	289	186	145	95	56	51	38
hbd_80 threshold	1914	1092	836	954	278	289	186	145	95	56	51	38
hbd_70 threshold	1933	1115	874	976	283	291	190	148	97	56	51	39
hbd_60 threshold	1928	1098	857	963	280	294	188	145	96	56	51	38
hbd_50 threshold	1933	1102	861	963	280	300	188	153	96	56	51	38
hbd_COG	1868	1117	860	777	282	291	151	148	97	55	49	42
thlA_90 threshold	1947	1216	1070	921	321	305	182	154	109	46	46	39
thlA_80 threshold	1971	1230	1086	927	324	305	184	154	109	48	48	39
thlA_70 threshold	2047	1248	1097	987	334	314	189	157	109	49	49	40
thlA_60 threshold	2036	1243	1094	981	334	313	187	159	109	49	49	40
thlA_50 threshold	2037	1263	1096	983	337	314	190	160	109	49	49	40
thlA_COG	1984	1155	1089	934	335	310	192	158	108	50	47	39", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Read 70 completeness total number data
total_data_raw <- read.table(text = "
Threshold	F.longum	F.prausnitzii	F.duncaniae	F.faecis	F.un1	F.un2	F.hattorii	F.un3	F.butyricigenerans	F.un4	F.wellingii	F.langellae
total_70	2614	1776	1448	1332	502	391	222	208	134	67	67	64", 
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define species
all_strains <- c("F.longum", "F.prausnitzii", "F.duncaniae", "F.faecis", "F.un1", "F.un2", 
                 "F.hattorii", "F.un3", "F.butyricigenerans", "F.un4", "F.wellingii", "F.langellae")

# Process gene data
gene_long <- gene_data_raw %>%
  separate(Gene, into = c("Gene", "Threshold"), sep = "_") %>%
  mutate(
    Gene = case_when(
      Gene == "cro" ~ "crt",
      TRUE ~ Gene
    ),
    Gene = factor(Gene, levels = c("thlA", "hbd", "crt", "bcd", "but")),
    Threshold = factor(Threshold, 
                       levels = c("90 threshold", "80 threshold", "70 threshold", 
                                  "60 threshold", "50 threshold", "COG"))
  ) %>%
  pivot_longer(
    cols = -c(Gene, Threshold),
    names_to = "Group",
    values_to = "Yes"
  ) %>%
  filter(Group %in% all_strains)

# Process total number data
total_long <- total_data_raw %>%
  pivot_longer(
    cols = -Threshold,
    names_to = "Group",
    values_to = "Total"
  ) %>%
  filter(Group %in% all_strains) %>%
  expand(nesting(Group, Total), 
         Threshold = c("90 threshold", "80 threshold", "70 threshold", 
                       "60 threshold", "50 threshold", "COG")) %>%
  mutate(
    Threshold = factor(Threshold, 
                       levels = c("90 threshold", "80 threshold", "70 threshold", 
                                  "60 threshold", "50 threshold", "COG"))
  )

# Merge data
all_data <- gene_long %>%
  left_join(total_long, by = c("Group", "Threshold")) %>%
  mutate(
    Yes = as.numeric(Yes),
    Total = as.numeric(Total),
    No = Total - Yes,
    Yes_percent = Yes / Total * 100,
    No_percent = No / Total * 100
  ) %>%
  select(Gene, Group, Threshold, Yes, No, Yes_percent, No_percent, Total)

# Convert data
plot_data_long <- all_data %>%
  pivot_longer(
    cols = c(Yes_percent, No_percent),
    names_to = "Category",
    values_to = "Percent"
  ) %>%
  mutate(
    Group = factor(Group, levels = all_strains),
    Category = ifelse(Category == "Yes_percent", "Yes", "No"),
    Count = ifelse(Category == "Yes", Yes, No),
    label = ifelse(Category == "Yes" & Percent > 0, sprintf("%.1f%%", Percent), ""),
    y_label = paste(Group, Threshold, sep = "_")
  ) %>%
  mutate(
    y_label = factor(y_label, levels = unlist(lapply(all_strains, function(strain) {
      paste(strain, c("90 threshold", "80 threshold", "70 threshold", 
                      "60 threshold", "50 threshold", "COG"), sep = "_")
    })))
  )

# color mapping
plot_data_long <- plot_data_long %>%
  mutate(
    fill_color = case_when(
      Category == "Yes" & Threshold == "90 threshold" ~ "#00441B",
      Category == "Yes" & Threshold == "80 threshold" ~ "#006D2C",
      Category == "Yes" & Threshold == "70 threshold" ~ "#238B45",
      Category == "Yes" & Threshold == "60 threshold" ~ "#41AB5D",
      Category == "Yes" & Threshold == "50 threshold" ~ "#74C476",
      Category == "Yes" & Threshold == "COG" ~ "#99CC33",
      Category == "No" & Threshold != "COG" ~ "#E5F5E0",
      Category == "No" & Threshold == "COG" ~ "#E5F5E0",
      TRUE ~ "#CCCCCC"
    )
  )

# Create plot by genes
create_gene_plot <- function(gene_name) {
  gene_data <- plot_data_long %>%
    filter(Gene == gene_name)

  n_strains <- length(all_strains)
  hline_positions <- seq(6.5, (n_strains-1)*6 + 0.5, by = 6)
  
  p <- ggplot(gene_data, aes(x = Percent, y = y_label, fill = fill_color)) +
    geom_bar(stat = "identity", width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = ifelse(Percent > 0.5, label, "")),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 2.0,  
              color = "white") +
    geom_hline(yintercept = hline_positions,
               linetype = "dashed", color = "black", size = 0.2) +
    scale_y_discrete(labels = NULL) +
    scale_x_continuous(
      limits = c(0, 100.01),
      breaks = seq(0, 100, 20),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    labs(title = gene_name, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey50", size = 0.5),
      plot.margin = margin(3, 3, 3, 3),
      legend.position = "none"
    ) +
    scale_fill_identity()
  
  return(p)
}
gene_names <- c("thlA", "hbd", "crt", "bcd", "but")
plots <- list()
for (gene in gene_names) {
  plots[[gene]] <- create_gene_plot(gene)
}

# Create species labels
create_strain_labels <- function() {
  n_strains <- length(all_strains)
  strain_positions <- seq(1.0, (n_strains-1)*6.5 + 1.0, by = 6.5)
  
  strain_labels <- data.frame(
    y = strain_positions,
    label = c(
      "italic(F.longum)", "italic(F.prausnitzii)", "italic(F.duncaniae)", 
      "italic(F.faecis)", "italic(F.un1)", "italic(F.un2)", 
      "italic(F.hattorii)", "italic(F.un3)", "italic(F.butyricigenerans)", 
      "italic(F.un4)", "italic(F.wellingii)", "italic(F.langellae)"
    )
  )
  
  max_y <- n_strains * 6 + 0.5
  
  p <- ggplot(strain_labels, aes(x = 0.5, y = y, label = label)) +
    geom_text(parse = TRUE, size = 2.5, fontface = "italic", hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, max_y), clip = "off")
  
  return(p)
}

strain_label_plot <- create_strain_labels()

# Create total number labels by species
create_total_labels_plot <- function() {
  n_strains <- length(all_strains)
  
  strain_totals <- total_data_raw %>%
    pivot_longer(
      cols = -Threshold,
      names_to = "Group",
      values_to = "Total"
    ) %>%
    filter(Group %in% all_strains) %>%
    mutate(
      Group = factor(Group, levels = all_strains)
    )
  
  strain_positions <- seq(1.0, (n_strains-1)*6.5 + 1.0, by = 6.5)
  
  total_labels <- data.frame(
    Group = all_strains,
    y_position = strain_positions,
    Total = strain_totals$Total[match(all_strains, strain_totals$Group)]
  )
  
  max_y <- n_strains * 6 + 0.5
  
  p <- ggplot(total_labels, aes(x = 0.5, y = y_position, 
                                label = paste0("n=", Total))) +
    geom_text(size = 2.5, hjust = 0.5, vjust = 0.5, fontface = "bold") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, max_y), clip = "off")
  
  return(p)
}

total_labels_plot <- create_total_labels_plot()

# Create legend
create_simple_legend <- function() {
  legend_color_mapping <- c(
    "90%" = "#00441B",
    "80%" = "#006D2C",
    "70%" = "#238B45",
    "60%" = "#41AB5D",
    "50%" = "#74C476",
    "COG" = "#99CC33"
  )
  
  legend_data <- data.frame(
    Threshold = names(legend_color_mapping),
    Color = legend_color_mapping,
    y = 6:1,
    x_color = 0.5,
    x_text = 0.7
  )
  
  p <- ggplot(legend_data) +
    geom_rect(aes(xmin = x_color - 0.15, xmax = x_color + 0.15, 
                  ymin = y - 0.35, ymax = y + 0.35, fill = Color),
              color = "black", size = 0.2) +
    geom_text(aes(x = x_text, y = y, label = Threshold), 
              vjust = 0.5, size = 2.5, color = "black", hjust = 0) +
    scale_fill_identity() +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold", 
                                margin = margin(b = 3))
    ) +
    labs(title = "Butyrate\nIdentity") +
    coord_cartesian(xlim = c(0.3, 1.2), ylim = c(0.5, 6.5), clip = "off")
  
  return(p)
}
simple_legend <- create_simple_legend()

create_empty_plot <- function() {
  p <- ggplot() + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  return(p)
}

empty_plot_top <- create_empty_plot()
empty_plot_bottom <- create_empty_plot()
legend_with_spaces <- wrap_plots(
  empty_plot_top,
  simple_legend,
  empty_plot_bottom,
  ncol = 1,
  heights = c(0.75, 0.2, 0.05)
)

# Combine plot
gene_plots <- (plots[["thlA"]] | plots[["hbd"]] | plots[["crt"]] | plots[["bcd"]] | plots[["but"]]) + 
  plot_layout(widths = rep(1, 5), nrow = 1)
final_plot <- wrap_plots(
  strain_label_plot,
  gene_plots,
  total_labels_plot,
  legend_with_spaces,
  ncol = 4,
  widths = c(1.2, 5, 0.4, 0.6)
) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

# Save as pdf
print(final_plot)
ggsave("butyrate_blasti_all.pdf", 
       plot = final_plot,
       width = 12,
       height = 12,
       units = "in",
       dpi = 300)
