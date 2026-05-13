# Import data
ani <- read.table("ani_matrix_typestrain_40.txt", sep = "\t")

# Load packages
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Convert and clean data
ani_t <- ani[, -c(4, 5)]
ani_long <- ani_t %>%
  spread(key = V1, value = V3)
row.names(ani_long) <- ani_long$V2
ani_long$V2 <- NULL

# Define typestrain
get_typestrain_name <- function(name) {
  if (startsWith(name, "A2165")) {
    return("F.duncaniae")
  } else if (startsWith(name, "GCA_036632915")) {
    return("F.taiwanense")
  } else if (startsWith(name, "APC922_41-1")) {
    return("F.hattorii")
  } else if (startsWith(name, "CLAJMH7B.combined")) {
    return("F.faecis")
  } else if (startsWith(name, "CLAAAH175.combined")) {
    return("F.tardum")
  } else if (startsWith(name, "CLAAAH281.combined")) {
    return("F.intestinale")
  } else if (startsWith(name, "CM04-06A")) {
    return("F.longum")
  } else if (startsWith(name, "AF52-21")) {
    return("F.butyricigenerans")
  } else if (startsWith(name, "ATCC_27768")) {
    return("F.prausnitzii")
  } else if (startsWith(name, "GCA_023347535")) {
    return("F.wellingii")
  } else if (startsWith(name, "GCF_002549775")) {
    return("F.langellae")
  } else if (startsWith(name, "GCA_041326985")) {
    return("F.un1")  
  } else if (startsWith(name, "GCA_025993375")) {
    return("F.un2")  
  } else if (startsWith(name, "GCA_041328755")) {
    return("F.un3")  
  } else if (startsWith(name, "GT_GMRSGB04542")) {
    return("F.un4")  
  } else if (startsWith(name, "GT_GMRSGB04546")) {
    return("F.un5")  
  } else if (startsWith(name, "GCA_934519045")) {
    return("F.un6")  
  } else if (startsWith(name, "GCA_959608795")) {
    return("F.un7")  
  } else if (startsWith(name, "GCA_934525985")) {
    return("F.un8")  
  } else if (startsWith(name, "MGYG000002651")) {
    return("F.un9")  
  } else if (startsWith(name, "GCA_934503275")) {
    return("F.un10")  
  } else if (startsWith(name, "GCA_934519035")) {
    return("F.un11")  
  } else if (startsWith(name, "GCA_934525735")) {
    return("F.un12")  
  } else if (startsWith(name, "GCA_934525845")) {
    return("F.un13")  
  } else if (startsWith(name, "GCA_934533135")) {
    return("F.un14")  
  } else if (startsWith(name, "GCA_938045145")) {
    return("F.un15")  
  } else if (startsWith(name, "GENOME000252")) {
    return("F.un16")  
  } else if (startsWith(name, "GCA_934501115")) {
    return("F.un17")  
  } else if (startsWith(name, "GCA_934553285")) {
    return("F.un18")  
  } else if (startsWith(name, "GCA_934526125")) {
    return("F.un19")  
  } else if (startsWith(name, "GCA_934548045")) {
    return("F.un20")  
  } else {
    return(NA)
  }
  #Faecalibacterium langellae：GCF_002549775.1
}

# Select typestrain
all_names <- union(rownames(ani_long), colnames(ani_long))
typestrain_flags <- !is.na(sapply(all_names, get_typestrain_name))
typestrain_names <- all_names[typestrain_flags]
ani_typestrain <- ani_long[
  rownames(ani_long) %in% typestrain_names,
  colnames(ani_long) %in% typestrain_names
]

# Set higher value
ani_typestrain_mat <- as.matrix(ani_typestrain)
ani_typestrain_mat <- pmax(ani_typestrain_mat, t(ani_typestrain_mat), na.rm = TRUE)

# Set typestrain name
name_mapping <- sapply(rownames(ani_typestrain_mat), get_typestrain_name)
rownames(ani_typestrain_mat) <- name_mapping
colnames(ani_typestrain_mat) <- name_mapping

# NA -> 80
ani_typestrain_mat[is.na(ani_typestrain_mat)] <- 80
diag(ani_typestrain_mat) <- 100

# Color mapping
heatmap_colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)

# Show 》95 value
ani_display <- round(ani_typestrain_mat, 1)
ani_display[ani_display < 95] <- NA
display_numbers_matrix <- ani_display
display_numbers_matrix[is.na(display_numbers_matrix)] <- ""

# Print and save
pdf("ANIheatmap_Typestrain_40.pdf", width = 8, height = 6)
pheatmap(ani_typestrain_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         fontsize_row = 6, 
         fontsize_col = 6,
         fontsize_number = 4,
         color = heatmap_colors,
         main = "ANI Matrix between Faecalibacterium typestrains",
         display_numbers = display_numbers_matrix,
         number_color = "black",
         cellwidth = 10, 
         cellheight = 10,
         border_color = NA,
         treeheight_row = 30,
         treeheight_col = 30
)
dev.off()

