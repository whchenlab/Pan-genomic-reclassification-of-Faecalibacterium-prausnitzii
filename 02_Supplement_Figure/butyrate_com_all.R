library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Import gene data
gene_data_raw <- read.table(text = "
Gene	Threshold	F.longum	F.prausnitzii	F.duncaniae	F.faecis	F.un1	F.un2	F.hattorii	F.un3	F.butyricigenerans	F.un4	F.wellingii	F.langellae
bcd_COG	com70	1765	1022	865	735	256	272	179	162	93	52	46	43
bcd_COG	com75	1698	993	825	712	224	257	176	129	88	51	43	42
bcd_COG	com80	1602	954	756	685	200	241	169	98	85	49	41	40
bcd_COG	com85	1416	890	695	650	166	220	162	80	80	45	38	39
bcd_COG	com90	1123	762	603	590	118	195	139	61	65	42	34	34
bcd_COG	com95	582	439	373	414	47	121	75	31	43	28	23	10
but_COG	com70	2414	1593	1331	1276	418	355	215	192	125	58	55	56
but_COG	com75	2301	1523	1245	1188	353	332	210	149	115	56	49	54
but_COG	com80	2124	1426	1130	1091	302	311	202	111	111	54	47	50
but_COG	com85	1841	1308	1023	984	232	289	192	87	104	49	45	49
but_COG	com90	1394	1084	860	850	166	252	161	66	84	44	39	42
but_COG	com95	670	542	466	520	57	137	83	34	47	29	24	16
cro_COG	com70	1934	1209	1040	878	263	323	179	153	104	54	50	39
cro_COG	com75	1855	1170	984	836	228	301	175	116	96	53	45	38
cro_COG	com80	1733	1108	900	792	200	283	168	90	93	51	43	38
cro_COG	com85	1531	1033	826	739	160	264	158	71	90	47	41	38
cro_COG	com90	1209	882	709	654	116	233	139	53	73	43	36	32
cro_COG	com95	610	490	419	439	46	131	75	29	42	29	24	14
hbd_COG	com70	1868	1117	860	777	282	291	151	148	97	55	49	42
hbd_COG	com75	1799	1081	820	748	250	269	150	115	90	54	45	40
hbd_COG	com80	1684	1032	747	718	222	253	144	88	87	52	43	38
hbd_COG	com85	1492	965	686	680	182	234	140	71	84	49	41	37
hbd_COG	com90	1180	827	599	612	133	204	120	53	67	45	36	31
hbd_COG	com95	601	473	370	420	53	128	66	29	40	29	23	11
thlA_COG	com70	1984	1155	1089	934	335	310	192	158	108	50	47	39
thlA_COG	com75	1901	1112	1029	888	287	287	188	121	99	49	43	38
thlA_COG	com80	1774	1061	934	843	251	269	180	94	96	47	41	38
thlA_COG	com85	1561	987	855	786	204	249	171	75	92	43	39	38
thlA_COG	com90	1227	833	731	701	149	219	149	57	74	40	34	32
thlA_COG	com95	618	471	420	465	56	122	80	30	43	28	24	15", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Import total number data
total_data_raw <- read.table(text = "
Threshold	F.longum	F.prausnitzii	F.duncaniae	F.faecis	F.un1	F.un2	F.hattorii	F.un3	F.butyricigenerans	F.un4	F.wellingii	F.langellae
total_95	680	556	479	523	59	140	83	35	49	30	27	16
total_90	1451	1150	904	858	190	267	166	68	91	49	45	44
total_85	1951	1410	1085	1004	271	311	198	90	113	55	53	55
total_80	2272	1559	1208	1119	360	339	208	116	120	61	58	58
total_75	2480	1678	1343	1226	420	363	217	159	124	63	61	62
total_70	2614	1776	1448	1332	502	391	222	208	134	67	67	64", 
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define species
all_strains <- c("F.longum", "F.prausnitzii", "F.duncaniae", "F.faecis", "F.un1", "F.un2", 
                 "F.hattorii", "F.un3", "F.butyricigenerans", "F.un4", "F.wellingii", "F.langellae")

# Process gene data
gene_long <- gene_data_raw %>%
  mutate(
    Gene = gsub("_COG", "", Gene),
    Gene = case_when(
      Gene == "cro" ~ "crt",
      TRUE ~ Gene
    ),
    Gene = factor(Gene, levels = c("thlA", "hbd", "crt", "bcd", "but")),
    Threshold = factor(Threshold, levels = c("com95", "com90", "com85", "com80", "com75", "com70"))
  ) %>%
  pivot_longer(
    cols = -c(Gene, Threshold),
    names_to = "Group",
    values_to = "Yes"
  ) %>%
  filter(Group %in% all_strains)

# Process total number data
total_long <- total_data_raw %>%
  mutate(
    Threshold = gsub("total_", "com", Threshold),
    Threshold = factor(Threshold, levels = c("com95", "com90", "com85", "com80", "com75", "com70"))
  ) %>%
  pivot_longer(
    cols = -Threshold,
    names_to = "Group",
    values_to = "Total"
  ) %>%
  filter(Group %in% all_strains)

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
      paste(strain, c("com95", "com90", "com85", "com80", "com75", "com70"), sep = "_")
    })))
  )

# Color mapping
color_mapping <- c(
  "com70_No" = "#EFEDF5",
  "com70_Yes" = "#BCBDDC",
  "com75_No" = "#EFEDF5",
  "com75_Yes" = "#9E9AC8",
  "com80_No" = "#EFEDF5",
  "com80_Yes" = "#807DBA",
  "com85_No" = "#EFEDF5",
  "com85_Yes" = "#6A51A3",
  "com90_No" = "#EFEDF5",
  "com90_Yes" = "#54278F",
  "com95_No" = "#EFEDF5",
  "com95_Yes" = "#3F007D"
)
plot_data_long <- plot_data_long %>%
  mutate(
    color_key = paste(Threshold, Category, sep = "_"),
    fill_color = color_mapping[color_key]
  )

# Create plots by genes
create_gene_plot <- function(gene_name) {
  gene_data <- plot_data_long %>%
    filter(Gene == gene_name)

  n_strains <- length(all_strains)
  hline_positions <- seq(6.5, (n_strains-1)*6 + 0.5, by = 6)
  
  p <- ggplot(gene_data, aes(x = Percent, y = y_label, fill = fill_color)) +
    geom_bar(stat = "identity", width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = ifelse(Percent > 1, label, "")),
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

# Create species number labels
create_total_labels_plot <- function() {
  n_strains <- length(all_strains)
  total_labels <- total_long %>%
    filter(Group %in% all_strains) %>%
    mutate(
      Group = factor(Group, levels = all_strains),
      strain_index = as.numeric(Group) - 1, 
      threshold_index = case_when(
        Threshold == "com95" ~ 0,
        Threshold == "com90" ~ 1,
        Threshold == "com85" ~ 2,
        Threshold == "com80" ~ 3,
        Threshold == "com75" ~ 4,
        Threshold == "com70" ~ 5
      ),
      y_position = strain_index * 6.6 + threshold_index - 2.2,
      label = paste0("n=", Total)
    ) %>%
    arrange(desc(strain_index), threshold_index)
  
  max_y <- n_strains * 6 + 0.5
  
  p <- ggplot(total_labels, aes(x = 0.5, y = y_position, label = label)) +
    geom_text(size = 1.8, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, max_y), clip = "off")
  
  return(p)
}

total_labels_plot <- create_total_labels_plot()

# Create legend
create_simple_legend <- function() {
  legend_color_mapping <- c(
    "com70" = "#BCBDDC",
    "com75" = "#9E9AC8",
    "com80" = "#807DBA",
    "com85" = "#6A51A3",
    "com90" = "#54278F",
    "com95" = "#3F007D"
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
    geom_text(aes(x = x_text, y = y, label = paste0(gsub("com", "", Threshold), "%")), 
              vjust = 0.5, size = 2.5, color = "black", hjust = 0) +
    scale_fill_identity() +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold", 
                                margin = margin(b = 3))
    ) +
    labs(title = "Butyrate\nCompleteness") +
    coord_cartesian(xlim = c(0.3, 1.2), ylim = c(0.5, 6.5), clip = "off")
  
  return(p)
}
simple_legend <- create_simple_legend()
legend_with_spaces <- wrap_plots(
  empty_plot_top,
  simple_legend,
  empty_plot_bottom,
  ncol = 1,
  heights = c(0.75, 0.2, 0.05)
)
create_empty_plot <- function() {
  p <- ggplot() + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  return(p)
}

# Combine plots
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
ggsave("butyrate_completeness_all.pdf", 
       plot = final_plot,
       width = 12,
       height = 12,
       units = "in",
       dpi = 300)
