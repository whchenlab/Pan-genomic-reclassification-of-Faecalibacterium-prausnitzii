# Import data
sgb <- read.table("sgb_matrix.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Load packages
library(tidyr)
library(pheatmap)
library(RColorBrewer)

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
}
sgb$Typestrain_Name <- sapply(sgb$V1, get_typestrain_name)
sgb <- sgb[!is.na(sgb$Typestrain_Name), ]

# Convert data
sgb_wide <- pivot_wider(
  sgb,
  id_cols = Typestrain_Name,
  names_from = V2,
  values_from = V3
)
sgb_matrix <- as.data.frame(sgb_wide)
rownames(sgb_matrix) <- sgb_matrix$Typestrain_Name
sgb_matrix$Typestrain_Name <- NULL
sgb_matrix_mat <- as.matrix(sgb_matrix)

# NA -> 85
sgb_matrix_mat[is.na(sgb_matrix_mat)] <- 85



# Color mapping
heatmap_colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
# Show >95 value
sgb_display <- round(sgb_matrix_mat, 1)
sgb_display[sgb_display < 95] <- NA
display_numbers_matrix <- sgb_display
display_numbers_matrix[is.na(display_numbers_matrix)] <- ""

# Save as pdf
pdf("ANIheatmap_SGB.pdf", width = 6, height = 6)
pheatmap(sgb_matrix_mat, 
         cluster_rows = TRUE,  
         cluster_cols = FALSE, 
         fontsize_row = 6, 
         fontsize_col = 6,
         fontsize_number = 4,
         color = heatmap_colors,
         main = " ",
         display_numbers = display_numbers_matrix,
         number_color = "black",
         cellwidth = 10, 
         cellheight = 10,
         border_color = NA,
         treeheight_row = 30,
         treeheight_col = 0,  
         show_colnames = TRUE,  
         show_rownames = TRUE   
)
dev.off()

