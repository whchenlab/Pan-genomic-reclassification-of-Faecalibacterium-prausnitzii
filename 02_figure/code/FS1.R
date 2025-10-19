determine_prefix <- function(names) {
  prefixes <- character(length(names))
  for (i in seq_along(names)) {
    if (startsWith(names[i], "A2165") || 
        startsWith(names[i], "GCA_032971305") || 
        startsWith(names[i], "GCA_010509575")) {
      prefixes[i] <- "F.duncaniae"
    } else if (startsWith(names[i], "GCA_036632915")) {
      prefixes[i] <- "F.taiwanense"
    } else if (startsWith(names[i], "APC922_41-1")) {
      prefixes[i] <- "F.hattorii"
    } else if (startsWith(names[i], "CLA-JM-H7B")) {
      prefixes[i] <- "F.faecis"
    } else if (startsWith(names[i], "CLA-AA-H175")) {
      prefixes[i] <- "F.tardum"
    } else if (startsWith(names[i], "CLA-AA-H281")) {
      prefixes[i] <- "F.intestinale"
    } else if (startsWith(names[i], "CM04-06A")) {
      prefixes[i] <- "F.longum"
    } else if (startsWith(names[i], "AF52-21")) {
      prefixes[i] <- "F.butyricigenerans"
    } else if (startsWith(names[i], "ATCC_27768")) {
      prefixes[i] <- "F.prausnitzii"
    } else if (startsWith(names[i], "F.wellingii_HTF-F")) {
      prefixes[i] <- "F.wellingii"
    } else if (startsWith(names[i], "F.gallinarum-JCM-17207")) {
      prefixes[i] <- "F.gallinarum"
    } else {
      prefixes[i] <- "Other"  # 其他情况
    }
  }
  return(prefixes)
}

# 确保相关列是字符型（避免因子型匹配问题）
Genome$Accession_Number <- as.character(Genome$Accession_Number)
Genome$Strain_Name <- as.character(Genome$Strain_Name)

# 分别基于Accession_Number和Strain_Name匹配类型
acc_types <- determine_prefix(Genome$Accession_Number)
strain_types <- determine_prefix(Genome$Strain_Name)

# 确定最终Type（优先用Accession_Number的结果，否则用Strain_Name的）
Genome$Type <- ifelse(acc_types != "Other", acc_types, strain_types)

# 提取符合条件的行，包含Accession_Number、Strain_Name、Source和Type列
Typestrain <- Genome[Genome$Type != "Other", c("Accession_Number", "Strain_Name", "Source", "Type")]

#保存为csv
write.csv(Typestrain, file = "Typestrain.csv", row.names = FALSE, fileEncoding = "UTF-8")

