# library rpackages ----
library('tidyverse')
library('ggsci')
library('ComplexHeatmap')
library('patchwork')
library('ggpubr')


# highlight.txt ----
highlight.txt = function(x,  color="red", family="") {
  #prevelence
  load('./01_data/oral_species_prevelence_PPI_baseline.RData')
  species_prevelence <- species_prevelence[which(species_prevelence$prevelence>0.1),]
  prevelence_s <- paste0('s__',rownames(species_prevelence))
  
  #ehomd
  eHOMD <- read.csv('./01_data/HOMD_taxon_table2023-05-03_1683123132.csv',header = T,skip = 1)
  oral_eHOMD <- subset(eHOMD,Body_site=='Oral')
  oral_eHOMD_s <- paste0('s__',oral_eHOMD$Genus,'_',oral_eHOMD$Species)
  oral_eHOMD_s <- unique(oral_eHOMD_s)
  
  
  
  
  #intersect_s
  intersect_s <- intersect(prevelence_s,oral_eHOMD_s)
  prevelence_s_diff <- setdiff(prevelence_s,intersect_s)
  oral_eHOMD_s_diff <- setdiff(oral_eHOMD_s,intersect_s)
  # library(glue)

  ha <- ifelse(x %in% intersect_s, glue("*#{x}"),
               ifelse(x %in% prevelence_s_diff, glue("*{x}"),
                      ifelse(x %in% oral_eHOMD_s_diff, glue("#{x}"),x)))
  return(ha)
}

highlight.txt.butyrate = function(x,  color="red", family="") {
  Butyrate <- get(load('01_data/Butyrate.RData'))
  Butyrate <- gsub(' ','_',Butyrate)
  Butyrate <- paste0('s__',Butyrate)
  # library(glue)

  x_t <- gsub('\\*|#','',x)
  
  ha <- ifelse(x_t %in% Butyrate, glue("~{x}"),x)
  return(ha)
}

highlight.txt.butyrate.new = function(x,  color="red", family="") {
  Butyrate <- get(load('01_data/Butyrate.new.RData'))
  Butyrate_genus <- subset(Butyrate,level=='genus')
  Butyrate_species <- subset(Butyrate,level=='species')
  
  Butyrate_genus <- Butyrate_genus$bac
  Butyrate_species <- Butyrate_species$bac
  # library(glue)
  
  
  x_t <- gsub('\\*|#','',x)
  

  
  for (bac in Butyrate_species) {
    x_t[grep(bac,x_t)] <- paste0('~',x_t[grep(bac,x_t)])
    
  }
  
  for (bac in Butyrate_genus) {
    x_t[grep(bac,x_t)] <- paste0('~',x_t[grep(bac,x_t)])
    
    # Butyrate_species <- Butyrate_species[-grep(bac,Butyrate_species)]
  }
  
  return(x_t)
}

highlight.txt.lactate.new = function(x,  color="red", family="") {
  Butyrate <- get(load('01_data/Lactate.new.RData'))
  Butyrate_genus <- subset(Butyrate,level=='genus')
  Butyrate_species <- subset(Butyrate,level=='species')
  
  Butyrate_genus <- Butyrate_genus$bac
  Butyrate_species <- Butyrate_species$bac
  # library(glue)
  
  
  x_t <- gsub('\\*|#','',x)
  
  for (bac in Butyrate_genus) {
    x_t[grep(bac,x_t)] <- paste0('!',x_t[grep(bac,x_t)])
    
    Butyrate_species <- Butyrate_species[-grep(bac,Butyrate_species)]
  }
  
  for (bac in Butyrate_species) {
    x_t[grep(bac,x_t)] <- paste0('!',x_t[grep(bac,x_t)])
    
  }
  
  return(x_t)
}

highlight.txt.propionate.new = function(x,  color="red", family="") {
  Butyrate <- get(load('01_data/Propionate.new.RData'))
  Butyrate_genus <- subset(Butyrate,level=='genus')
  Butyrate_species <- subset(Butyrate,level=='species')
  
  Butyrate_genus <- Butyrate_genus$bac
  Butyrate_species <- Butyrate_species$bac
  # library(glue)
  
  
  x_t <- gsub('\\*|#','',x)
  
  for (bac in Butyrate_genus) {
    x_t[grep(bac,x_t)] <- paste0('+',x_t[grep(bac,x_t)])
    
    Butyrate_species <- Butyrate_species[-grep(bac,Butyrate_species)]
  }
  
  for (bac in Butyrate_species) {
    x_t[grep(bac,x_t)] <- paste0('+',x_t[grep(bac,x_t)])
    
  }
  
  return(x_t)
}

highlight.txt.acetate.new = function(x,  color="red", family="") {
  Butyrate <- get(load('01_data/Acetate.new.RData'))
  Butyrate_genus <- subset(Butyrate,level=='genus')
  Butyrate_species <- subset(Butyrate,level=='species')
  
  Butyrate_genus <- Butyrate_genus$bac
  Butyrate_species <- Butyrate_species$bac
  # library(glue)
  
  
  x_t <- gsub('\\*|#','',x)
  
  for (bac in Butyrate_genus) {
    x_t[grep(bac,x_t)] <- paste0('-',x_t[grep(bac,x_t)])
    
    Butyrate_species <- Butyrate_species[-grep(bac,Butyrate_species)]
  }
  
  for (bac in Butyrate_species) {
    x_t[grep(bac,x_t)] <- paste0('-',x_t[grep(bac,x_t)])
    
  }
  
  return(x_t)
}

# namediff ----
namediff <- function(DF0,DF1){
  # DF0 <- feat_t
  colname0 <- colnames(DF0)
  # DF1 <- feat_unique
  colname1 <- colnames(DF1)
  diffname <- setdiff(colname0,colname1)
  # grep('^s__St',a,value = T)
  # grep('^s__St',colname0,value = T)
  # grep('^s__St',colname1,value = T)
  DF1_t <- as.data.frame(matrix(data=0,nrow = nrow(DF1),ncol = length(diffname)))
  colnames(DF1_t) <- diffname
  DF1 <- cbind(DF1,DF1_t)
  DF1 <- DF1[,colname0]
  return(DF1)
}
# metaphlan4.species ----
#row is species
metaphlan4.species <- function(data,name_exclude=T){
  data_species <- data[-grep('t__',data$clade_name,value=FALSE),]
  data_species <- data_species[grep('s__',data_species$clade_name,value=FALSE),]
  data_species$clade_name <- sapply(strsplit(data_species$clade_name,'\\|'),'[',7)
  data_species$clade_name <- sub("s__", "", data_species$clade_name)
  data_species <- as.data.frame(t(data_species))
  colnames(data_species) <- data_species['clade_name',]
  data_species <- data_species[-1,]

  if(length(unique(colnames(data_species)))<ncol(data_species)){
    message('dumplicated species!')
    species_unique <- unique(colnames(data_species))
    feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(data_species),ncol = length(species_unique)))
    colnames(feat_unique) <- species_unique
    rownames(feat_unique) <- rownames(data_species)
    for (species_t in species_unique) {
      match_t <- colnames(data_species)==species_t
      if(sum(match_t)>1){
        feat_t <- data_species[,colnames(data_species)==species_t]
        feat_unique[,species_t] <- rowSums(feat_t)
      }else{
        feat_unique[,species_t] <- data_species[,species_t]
      }
    }
    data_species <- feat_unique
  }
  # rowSums(feat_unique)
  # rowSums(feat)
  if(name_exclude==T){
  rownames(data_species) <- sapply(strsplit(rownames(data_species),'\\.'),'[',1)
  }
  data_species <- as.data.frame(t(data_species))
  # library(dplyr)
  data_species %>% mutate(across(where(is.character), as.numeric)) -> data_species
  return(data_species)
}

# mOTU.species ----
#row is species
mOTU.species <- function(data,name_exclude=T){
  # data_species <- data[-grep('t__',data$clade_name,value=FALSE),]
  data_species <- data
  data_species <- data_species[grep('s__',data_species$clade_name,value=FALSE),]
  data_species$clade_name <- sapply(strsplit(data_species$clade_name,'\\|'),'[',7)
  data_species$clade_name <- sub("s__", "", data_species$clade_name)
  
  data_species <- as.data.frame(t(data_species),check.names = FALSE)
  colnames(data_species) <- data_species['clade_name',]
  data_species <- data_species[-1,]
  
  if(length(unique(colnames(data_species)))<ncol(data_species)){
    message('dumplicated species!')
    species_unique <- unique(colnames(data_species))
    feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(data_species),ncol = length(species_unique)))
    colnames(feat_unique) <- species_unique
    rownames(feat_unique) <- rownames(data_species)
    for (species_t in species_unique) {
      match_t <- colnames(data_species)==species_t
      if(sum(match_t)>1){
        feat_t <- data_species[,colnames(data_species)==species_t]
        feat_t %>% mutate(across(where(is.character), as.numeric)) -> feat_t
        feat_unique[,species_t] <- rowSums(feat_t)
      }else{
        feat_unique[,species_t] <- data_species[,species_t]
      }
    }
    data_species <- feat_unique
  }
  # rowSums(feat_unique)
  # rowSums(feat)
  # if(name_exclude==T){
  #   rownames(data_species) <- sapply(strsplit(rownames(data_species),'\\.'),'[',1)
  # }
  data_species <- t(data_species)
  data_species <- as.data.frame(data_species,check.names = FALSE)
  # library(dplyr)
  data_species %>% mutate(across(where(is.character), as.numeric)) -> data_species
  return(data_species)
}

# curatedMetagenomicData.species ----
#row is species
curatedMetagenomicData.species <- function(data,name_exclude=T){
  data_species <- as.data.frame(t(data))
  data_species <- data_species[grep('s__', rownames(data_species),value=FALSE),]
  data_species$clade_name <- sapply(strsplit(rownames(data_species),'\\.'),'[',7)
  data_species$clade_name <- sub("s__", "", data_species$clade_name)
  data_species <- as.data.frame(t(data_species))
  colnames(data_species) <- data_species['clade_name',]
  data_species <- data_species[-which(rownames(data_species)=='clade_name'),]
  
  if(length(unique(colnames(data_species)))<ncol(data_species)){
    message('dumplicated species!')
    species_unique <- unique(colnames(data_species))
    feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(data_species),ncol = length(species_unique)))
    colnames(feat_unique) <- species_unique
    rownames(feat_unique) <- rownames(data_species)
    for (species_t in species_unique) {
      match_t <- colnames(data_species)==species_t
      if(sum(match_t)>1){
        feat_t <- data_species[,colnames(data_species)==species_t]
        feat_unique[,species_t] <- rowSums(feat_t)
      }else{
        feat_unique[,species_t] <- data_species[,species_t]
      }
    }
    data_species <- feat_unique
  }
  # rowSums(feat_unique)
  # rowSums(feat)
  if(name_exclude==T){
    rownames(data_species) <- sapply(strsplit(rownames(data_species),'\\.'),'[',1)
  }
  data_species <- as.data.frame(t(data_species))
  # library(dplyr)
  data_species %>% mutate(across(where(is.character), as.numeric)) -> data_species
  return(data_species)
}

# metaphlan4.genus ----
metaphlan4.genus <- function(data,name_exclude=T){
  data_species <- data[-grep('t__',data$clade_name,value=FALSE),]
  data_species <- data_species[grep('g__',data_species$clade_name,value=FALSE),]
  data_species$clade_name <- sapply(strsplit(data_species$clade_name,'\\|'),'[',6)
  # data_species$clade_name <- sub("g__", "", data_species$clade_name)
  data_species <- as.data.frame(t(data_species))
  colnames(data_species) <- data_species['clade_name',]
  data_species <- data_species[-1,]
  
  if(length(unique(colnames(data_species)))<ncol(data_species)){
    message('dumplicated species!')
    species_unique <- unique(colnames(data_species))
    feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(data_species),ncol = length(species_unique)))
    colnames(feat_unique) <- species_unique
    rownames(feat_unique) <- rownames(data_species)
    for (species_t in species_unique) {
      match_t <- colnames(data_species)==species_t
      if(sum(match_t)>1){
        feat_t <- data_species[,colnames(data_species)==species_t]
        feat_t <- apply(feat_t, 2, as.numeric)
        feat_unique[,species_t] <- rowSums(feat_t)
      }else{
        feat_t <- apply(feat_t, 2, as.numeric)
        feat_unique[,species_t] <- data_species[,species_t]
      }
    }
    data_species <- feat_unique
  }
  # rowSums(feat_unique)
  # rowSums(feat)
  if(name_exclude==T){
    rownames(data_species) <- sapply(strsplit(rownames(data_species),'\\.'),'[',1)
  }
  data_species <- as.data.frame(t(data_species))
  # library(dplyr)
  data_species %>% mutate(across(where(is.character), as.numeric)) -> data_species
  return(data_species)
}

# abundance_norm ----
# col is sample
#' @param data col is sample, row is feature
#' @param nor if last possess is normalized
#' @param filter.taxon filter taxon or not
#' @param percent if last possess is percent
#' @param filter.sample filter sample or not
#' @param metaphlan from metaphlan or not
abundance.norm <- function(data,nor=T,filter.taxon=T,percent=F,filter.sample=F,metaphlan=F){
  genus.data <- data
  if(metaphlan==F){
    genus.data <- genus.data/matrix(colSums(genus.data), nrow = nrow(genus.data) , ncol = ncol(genus.data) , byrow=TRUE)  
  }else{
    genus.data <- genus.data/100    
  }
  if(filter.taxon==T){
    ori_len <- nrow(genus.data)
    #exclude prevelence < 2
    genus.data <- genus.data[rowSums(genus.data!=0)>2,]
    #exclude all abundance < 0.001
    genus.data <- genus.data[rowSums(genus.data>=0.001)>0,]
    now_len <- nrow(genus.data)
    print(paste0('++++ filter ',ori_len-now_len,' taxon!'))
    print(paste0('++++ retain ',now_len,' taxon!'))
    
  }
  if(filter.sample==T){
    #exclude prevelence < 2
    ori_len <- ncol(genus.data)
    genus.data <- genus.data[,colSums(genus.data!=0)>2]
    now_len <- ncol(genus.data)
    print(paste0('++++ filter ',ori_len-now_len,' samples!'))
    print(paste0('++++ retain ',now_len,' samples!'))
  }
  #normalization
  if(nor==T){
    genus.data <- genus.data/matrix(colSums(genus.data), nrow = nrow(genus.data) , ncol = ncol(genus.data) , byrow=TRUE)
  }
  if(percent==T){
    genus.data <- genus.data*100
  }
  genus.data <- as.data.frame(t(genus.data))
  return(genus.data)
}



# R for stamp ----
#data is a dataframe: col is sample, row is feature; group is a dataframe, which is paided with data.
#alternative is two.sided or one
Rstamp <- function(data,group,method='wilcox',p_adjust=T,p.matrix=F,alternative="two.sided",only.matriax=F, plot.sort="p.value",
                   filter0=1,paired=F,confidence=T){
  
  # library(tidyverse)
  # data <- data*100
  # data <- data %>% filter(apply(data,1,mean) > 1)
  data <- t(data)
  data1 <- data.frame(data,group$group)
  colnames(data1) <- c(colnames(data),"Group")
  data1$Group <- as.factor(data1$Group)
  if(method=='wilcox'){
    
    if(alternative=="two.sided"){
      if(paired==F){
        diff <- data1 %>%
        select_if(is.numeric) %>%
        map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1,conf.int = T)), .id = 'var')
      }else{
        diff <- data1 %>%
          select_if(is.numeric) %>%
          map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1,conf.int = T,paired = T)), .id = 'var')
      }
      
    }else{
      diff <- data.frame()
      diff1 <- data1 %>%
        select_if(is.numeric) %>%
        map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1,conf.int = T,alternative = "less")), .id = 'var')
      diff2 <- data1 %>%
        select_if(is.numeric) %>%
        map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1,conf.int = T,alternative = "greater")), .id = 'var')
      for (diff_i in 1:nrow(diff1)) {
        p1 <- diff1$p.value[diff_i]
        p2 <- diff2$p.value[diff_i]
        if(p1<p2){
          diff <- rbind(diff,diff1[diff_i,])
        }else{
          diff <- rbind(diff,diff2[diff_i,])
        }
      }
    }
    # }else if(alternative=="less"){
    #   diff <- data1 %>%
    #     select_if(is.numeric) %>%
    #     map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1,conf.int = T,alternative = "less")), .id = 'var')
    # }else if(alternative=="greater"){
    #   diff <- data1 %>%
    #     select_if(is.numeric) %>%
    #     map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data1,conf.int = T,alternative = "greater")), .id = 'var')
    # }
    
  }else{
    diff <- data1 %>%
      select_if(is.numeric) %>%
      map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')
  }
  
  if(p_adjust==T){
    diff$p.value <- p.adjust(diff$p.value,"BH")
  }
  
  
  
  if(only.matriax==T){
    if(plot.sort=="estimate"){
      diff=diff[sort(diff$estimate,index.return = T)$ix,]
    }else{
      diff=diff[sort(diff$p.value,index.return = T)$ix,]
    }
    return(diff)
  }
  
  diff <- diff %>% filter(p.value <= filter0)
  
  abun.bar <- data1[,c(diff$var,"Group")] %>%
    gather(variable,value,-Group) %>%
    group_by(variable,Group) %>%
    summarise(Mean = mean(value))
  diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                              levels(data1$Group)[2]))
  if(plot.sort=="estimate"){
    diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]
  }else{
    diff.mean <- diff.mean[order(diff.mean$p.value,decreasing = F),]
  }
  
  
  #plot
  # library(ggplot2)
  cbbPalette <- c("#91B7AC","#EFBC9A")
  abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))
  
  p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
    scale_x_discrete(limits = levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),#margin = margin(r = 20)
          legend.title=element_blank(),
          legend.text=element_text(size=12,face = "bold",colour = "black"),
          legend.position = c(-1,-0.1),
          legend.direction = "horizontal",
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.5,"cm"))
  
  
  for (i in 1:(nrow(diff.mean) - 1)){
    p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }
  
  p1 <- p1 +
    geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
    scale_fill_manual(values=cbbPalette)
  # p1
  
  
  
  #
  
  diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
  diff.mean$p.value <- signif(diff.mean$p.value,3)
  diff.mean$p.value <- as.character(diff.mean$p.value)
  
  p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
    scale_x_discrete(limits = levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals")
  
  for (i in 1:(nrow(diff.mean) - 1)){
    p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }
  
  p2 <- p2 +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.8), width = 0.5, size = 0.5) +
    geom_point(shape = 21,size = 3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')
  
  # p2
  
  #
  if(p_adjust==T){
    p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
      geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
      geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
                srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }else{
    p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
      geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
                hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
      geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value",
                srt = 90,fontface = "bold",size = 5) +
      coord_flip() +
      ylim(c(0,1)) +
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())
  }
  # p3
  
  # library(patchwork)
  if(confidence==T){
    p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))
  }else{
    p <- p1 + p3 + plot_layout(widths = c(4,2))
  }
  
  
  if(p.matrix==T){
    Rstamp.result <- list(p.matrix=diff,p=p)
    return(Rstamp.result)
  }else{
    return(p)
  }
}

# qiime2 ----
qiime2.genus <- function(genus.data0,nor=F){
  genus.data <- genus.data0[-grep('_$|uncultured$|Incertae_Sedis$|(gut metagenome$)|(uncultured bacterium$)|(uncultured organism$)',rownames(genus.data0),value=FALSE),]
  if(sum(table(paste0(sapply(strsplit(rownames(genus.data),';'),"[",6)))>1)!=0){
    print('There are dumplicated genus!')
    Strgenus <- names(which(table(paste0(sapply(strsplit(rownames(genus.data),';'),"[",6)))>1))
    feat1 <- genus.data[!sapply(strsplit(rownames(genus.data),';'),"[",6) %in% Strgenus,]
    feat3 <- data.frame()
    for (i in Strgenus) {
      feat2 <- genus.data[sapply(strsplit(rownames(genus.data),';'),"[",6) %in% i,]
      feat2 <- t(data.frame(colSums(feat2)))
      row.names(feat2) <- i
      feat3 <- rbind(feat3,feat2)
    }
    rownames(feat1) <- paste0(sapply(strsplit(rownames(feat1),';'),"[",6))
    genus.data <- rbind(feat1,feat3)
  }else{
    rownames(genus.data) <- paste0(sapply(strsplit(rownames(genus.data),';'),"[",6))
  }
  #exclude prevelence < 2
  genus.data <- genus.data[rowSums(genus.data!=0)>2,]
  #normalization
  if(nor==T){
    genus.data <- genus.data/matrix(colSums(genus.data), nrow = nrow(genus.data) , ncol = ncol(genus.data) , byrow=TRUE)
  }
  
  genus.data <- t(genus.data)
  return(genus.data)
}
# marker plot ----
marker_plot_heatmap <- function(marker1,
                                highlight=F,
                                bar_griup=F,
                                cluster_rows_is=F,
                                cluster_columns_is=F){
  # library(tidyr)
  # library(ggplot2)
  # library(ComplexHeatmap)
  # library(RColorBrewer)
  
  
  marker1 <- marker1[,c("scientific_name","project_id","LDA")]
  marker1$LDA <- as.numeric(marker1$LDA)
  plot.data.1 <- marker1 %>% spread(project_id,LDA)
  rownames(plot.data.1) <- plot.data.1$scientific_name
  plot.data.1$scientific_name <- NULL
  plot.data.1[is.na(plot.data.1)] <- 0
  # rname <- colnames(plot.data.1)
  # colnames(plot.data.1) <- paste0('C',1:ncol(plot.data.1))
  # plot.data.1 <- arrange(plot.data.1,C2,C1)
  # colnames(plot.data.1) <- rname
  plot.data.1[plot.data.1==0] <- NA
  
  plot.data.2 <- plot.data.1
  plot.data.2[is.na(plot.data.2)] <- 0
  plot.data.2 <- data.frame('Disease' = rowSums(plot.data.2>0),'Health' = rowSums(plot.data.2<0))
  
  plot.data.1 <- cbind(plot.data.1,plot.data.2)
  plot.data.1 <- arrange(plot.data.1,Disease,desc(Health))
  C1 <- ncol(plot.data.1)
  C2 <- ncol(plot.data.1)-1
  plot.data.1 <- plot.data.1[,-c(C1,C2)]
  
  if(bar_griup==F){
    bar1 <- HeatmapAnnotation(
      sum1 = anno_barplot(
        colSums(!is.na(plot.data.1)),
        bar_width = 0.9,
        gp = gpar(col = "white", fill = "#8988A3"),
        border = F,
        # axis_param = list(at = c(0,1.25e5,2.5e5),
        #                   labels = c("0","125k","250k")),
        height = unit(2, "cm")
      ), show_annotation_name = F
    )
    
    bar2 <- rowAnnotation(
      sum2 = anno_barplot(
        rowSums(!is.na(plot.data.1)),
        bar_width = 0.9,
        gp = gpar(col = "white", fill = "#8988A3"),
        border = F,
        # axis_param = list(at = c(0,2.5e5,5e5),
        #                   labels = c("0","250k","500k")),
        width = unit(2, "cm")
      ), show_annotation_name = F
    )}else{
      fill_cl <- c('#cf9198','#c4c7b4')#'#cf9198','#c4c7b4'  '#C41E26','#313695'  '#D13737','#6298C5'  '#ecd6bc','#c4c7b4'
      plot.data.2 <- plot.data.1
      plot.data.2[is.na(plot.data.2)] <- 0
      plot.data.2 <- data.frame('Disease' = colSums(plot.data.2>0),'Health' = colSums(plot.data.2<0))
      bar1 <- HeatmapAnnotation(
        sum1 = anno_barplot(
          plot.data.2,
          bar_width = 0.9,
          gp = gpar(col = "white", fill = fill_cl),
          border = F,
          # axis_param = list(at = c(0,1.25e5,2.5e5),
          #                   labels = c("0","125k","250k")),
          height = unit(2, "cm")
        ), show_annotation_name = F
      )
      # draw(bar1)
      
      plot.data.2 <- plot.data.1
      plot.data.2[is.na(plot.data.2)] <- 0
      plot.data.2 <- data.frame('Disease' = rowSums(plot.data.2>0),'Health' = rowSums(plot.data.2<0))
      bar2 <- rowAnnotation(
        sum2 = anno_barplot(
          plot.data.2,
          bar_width = 0.9,
          gp = gpar(col = "white", fill = fill_cl), ##aa3e53 #48597e
          border = F,
          # axis_param = list(at = c(0,2.5e5,5e5),
          #                   labels = c("0","250k","500k")),
          # width = unit(2, "cm")
        ), show_annotation_name = F
      )
      # draw(bar2)
    }
  library(circlize)
  col_fun = colorRamp2(c(-3,0,3),c("#313695", "white", "#A50026"))
  # col_fun = colorRamp2(c(-3,0,3),c("#597C8B", "white", '#b36a6f'))
  # col_fun(seq(-3, 3))
   
  if(highlight==T){
    p <- Heatmap(as.matrix(plot.data.1),
                 # col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
                 col = col_fun,
                 na_col = "white",
                 row_labels = highlight.txt.butyrate(highlight.txt(rownames(plot.data.1))),
                 # row_title=1,
                 row_names_side = "left",
                 cluster_rows = cluster_rows_is,
                 cluster_columns = cluster_columns_is,
                 rect_gp = gpar(col = "white", lwd = 2),
                 name = "~ Butyrate\n\n* Oral self\n\n# ehomd\n\nLDA score",
                 top_annotation = bar1,
                 right_annotation = bar2,
                 # left_annotation = bar3,
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 width = ncol(plot.data.1) * unit(2.5, "mm"),
                 height = nrow(plot.data.1) * unit(2.5, "mm"))
    
  }else{
    p <- Heatmap(as.matrix(plot.data.1),
                 # col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
                 col = col_fun,
                 na_col = "white",
                 row_names_side = "left",
                 cluster_rows = cluster_rows_is,
                 cluster_columns = cluster_columns_is,
                 rect_gp = gpar(col = "white", lwd = 2),
                 name = "LDA score",
                 top_annotation = bar1,
                 right_annotation = bar2,
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 width = ncol(plot.data.1) * unit(2.5, "mm"),
                 height = nrow(plot.data.1) * unit(2.5, "mm")
    )
  }
  
  
  p
  
}

# calculate.FPR ----
# 1 is the cutoff was selected by train model, and 2 is the cutoff was selected by testing data.
#' @param model_siamcat is model$models_top$PPI_self, model is the output of my_result or my_result_adj
#' @param self.ext is model$cross_models$`PPI_self-MetaCarkis_3`
calculate.FPR <- function(model_siamcat,self.ext,return_label=F,threshold_standard='F1_score'){
  # 1
  siamcat.holdout <- self.ext$model_back
  models_Y <- model_siamcat$model_back
  pred <- rowMeans(pred_matrix(siamcat.holdout))
  if(threshold_standard=='F1_score'){
    sensitivities <- models_Y@eval_data$roc$sensitivities
    specificities <- models_Y@eval_data$roc$specificities
    F1_scores <- 2*(sensitivities*specificities)/(sensitivities+specificities)
    threshold1 <- models_Y@eval_data$roc$thresholds[which(F1_scores==max(F1_scores))[1]]
  }else{
    threshold1 <- models_Y@eval_data$roc$thresholds[which(models_Y@eval_data$roc$specificities > 0.90)[1]]
  }
  
  df.fpr <- tibble(pred=pred, Group=(label(siamcat.holdout))[[1]])
  df.fpr.df1 <- data.frame(pred=pred, Group=(label(siamcat.holdout))[[1]],row.names = names(pred))
  df.fpr.df1$Label <- 'Train'
  df.fpr.df1$threshold <- threshold1
  df.fpr.df1$pred_label <- ifelse(df.fpr.df1$pred>threshold1,1,-1)
  df.fpr.df1$sample <- rownames(df.fpr.df1)
  
  df.fpr.temp1 <- df.fpr %>% 
    group_by(Group) %>% 
    dplyr::summarise(n=n(), 
                     pred.pos=sum(pred>threshold1), 
                     pred.neg=sum(pred < threshold1)) %>% 
    ungroup()
  train_TPR <- df.fpr.temp1[which(df.fpr.temp1$Group==1),'pred.pos']/df.fpr.temp1[which(df.fpr.temp1$Group==1),'n']
  train_FPR <- df.fpr.temp1[which(df.fpr.temp1$Group==-1),'pred.pos']/df.fpr.temp1[which(df.fpr.temp1$Group==-1),'n']
  # knitr::kable(df.fpr.temp1)
  train_mat <- df.fpr.temp1
  confonder_frame <- cbind(train_mat,Label=c("FPR","TPR"),num=c(train_FPR[1,1],train_TPR[1,1]),Class='Train')
  
  # 2
  siamcat.holdout <- self.ext$model_back
  models_Y <- self.ext$model_back
  pred <- rowMeans(pred_matrix(siamcat.holdout))
  if(threshold_standard=='F1_score'){
    sensitivities <- models_Y@eval_data$roc$sensitivities
    specificities <- models_Y@eval_data$roc$specificities
    F1_scores <- 2*(sensitivities*specificities)/(sensitivities+specificities)
    threshold2 <- models_Y@eval_data$roc$thresholds[which(F1_scores==max(F1_scores))[1]]
  }else{
    threshold2 <- models_Y@eval_data$roc$thresholds[which(models_Y@eval_data$roc$specificities > 0.90)[1]]
  }
  
  df.fpr <- tibble(pred=pred, Group=(label(siamcat.holdout))[[1]])
  df.fpr.df2 <- data.frame(pred=pred, Group=(label(siamcat.holdout))[[1]],row.names = names(pred))
  df.fpr.df2$Label <- 'Test'
  df.fpr.df2$threshold <- threshold2
  df.fpr.df2$pred_label <- ifelse(df.fpr.df1$pred>threshold2,1,-1)
  df.fpr.df2$sample <- rownames(df.fpr.df2)
  
  df.fpr.temp1 <- df.fpr %>% 
    group_by(Group) %>% 
    dplyr::summarise(n=n(), 
                     pred.pos=sum(pred>threshold2), 
                     pred.neg=sum(pred < threshold2)) %>% 
    ungroup()
  test_TPR <- df.fpr.temp1[which(df.fpr.temp1$Group==1),'pred.pos']/df.fpr.temp1[which(df.fpr.temp1$Group==1),'n']
  test_FPR <- df.fpr.temp1[which(df.fpr.temp1$Group==-1),'pred.pos']/df.fpr.temp1[which(df.fpr.temp1$Group==-1),'n']
  # knitr::kable(df.fpr.temp2)
  test_mat <- df.fpr.temp1
  confonder_frame <- rbind(confonder_frame,
                           cbind(test_mat,Label=c("FPR","TPR"),num=c(test_FPR[1,1],test_TPR[1,1]),Class='Test'))
  
  
  # 3
  siamcat.holdout <- self.ext$model_back
  pred <- rowMeans(pred_matrix(siamcat.holdout))
 
  threshold3 <- mean(c(threshold1,threshold2))
  
  df.fpr <- tibble(pred=pred, Group=(label(siamcat.holdout))[[1]])
  df.fpr.df3 <- data.frame(pred=pred, Group=(label(siamcat.holdout))[[1]],row.names = names(pred))
  df.fpr.df3$Label <- 'Mean'
  df.fpr.df3$threshold <- threshold3
  df.fpr.df3$pred_label <- ifelse(df.fpr.df1$pred>threshold3,1,-1)
  df.fpr.df3$sample <- rownames(df.fpr.df3)
  
  df.fpr.temp1 <- df.fpr %>% 
    group_by(Group) %>% 
    dplyr::summarise(n=n(), 
                     pred.pos=sum(pred>threshold3), 
                     pred.neg=sum(pred < threshold3)) %>% 
    ungroup()
  test_TPR <- df.fpr.temp1[which(df.fpr.temp1$Group==1),'pred.pos']/df.fpr.temp1[which(df.fpr.temp1$Group==1),'n']
  test_FPR <- df.fpr.temp1[which(df.fpr.temp1$Group==-1),'pred.pos']/df.fpr.temp1[which(df.fpr.temp1$Group==-1),'n']
  # knitr::kable(df.fpr.temp2)
  test_mat <- df.fpr.temp1
  confonder_frame <- rbind(confonder_frame,
                           cbind(test_mat,Label=c("FPR","TPR"),num=c(test_FPR[1,1],test_TPR[1,1]),Class='Mean'))
  
  
  pred_df <- rbind(df.fpr.df1,df.fpr.df2,df.fpr.df3)
  if(return_label==T){
    return(pred_df)
  }else{
    return(confonder_frame)
  }
}




# plot_FPR_TPR ----
plot_FPR_TPR <- function(confonder_frame_all,y_label0,figure_name,figure_number='',return_label=F){
  # plot.bar
  # library(ggplot2)
  # library(ggpubr)
  # library(cowplot)
  # library(RColorBrewer)
  # library(dplyr)
  # library(tidyr)
  # library(xlsx2dfs)
  # library(ggsignif)
  

  if(return_label==T){
    confonder_frame_all <- subset(confonder_frame_all,Label==figure_name)
    
    if(y_label0=='TPR'){
      Group0=1
      figure_name <- paste0(figure_name,'; ','Disease')
    }else{
      Group0=-1
      figure_name <- paste0(figure_name,'; ','Health')
    }
    confounder_bar_PPI <- subset(confonder_frame_all,Group==Group0) %>% 
      group_by(label) %>%  
      summarise(PPI=sum(pred_label==1),
                Non_PPI=sum(pred_label==-1))
    data_plot <- confounder_bar_PPI %>% gather('Group','Num',-label)
    
    p <- ggplot(data_plot, aes(label, weight = Num)) +
      # geom_hline(yintercept = seq(10, 50, 10), color = 'gray') +
      geom_bar(aes(fill = Group),color = "white", width = .65, position = 'dodge') +
      # scale_fill_manual(palette ="Set3")+  #values = c('#90BFF9','#05BE78','#FF0000','#B07FC0','#C5944E')
      scale_fill_manual(values = c('#7DABCF','#F46149'))+ #AAB083 939650 4D7787 7DABCF   D6594C
      # geom_errorbar(aes(ymin = num - se, ymax = num + se), width = 0.25, size = 0.3, position = position_dodge(0.7)) +
      labs( y = 'No. of sample',x='') +
      # ylim(0,1)+
      # scale_y_continuous(expand = c(0,0),limits = c(0, length0*1.55)) +
      geom_hline(yintercept =10,color='#EE4431')+
      geom_hline(yintercept =5,color='blue')+
      ggtitle(figure_number,figure_name)+
      theme_classic()+
      theme(legend.position="right",
            axis.text.x=element_text(angle=45, hjust=1,face = 'plain',size=10),
            axis.text.y=element_text(face = 'plain',size=10),
            text = element_text(size=15,face = 'plain',family ="sans",colour = 'black'))+
      theme(plot.title = element_text(hjust = -0.1, vjust=1,size=18))
  }else{
    data_plot <- subset(confonder_frame_all,(Class==figure_name)&(Label==y_label0))
    if(y_label0=='FPR'){
      p <- ggplot(data_plot, aes(label, weight = 1-num, fill = label)) +
        # geom_hline(yintercept = seq(10, 50, 10), color = 'gray') +
        geom_bar(color = "white", width = .65, position = 'dodge') +
        # scale_fill_manual(palette ="Set3")+  #values = c('#90BFF9','#05BE78','#FF0000','#B07FC0','#C5944E')
        # geom_errorbar(aes(ymin = num - se, ymax = num + se), width = 0.25, size = 0.3, position = position_dodge(0.7)) +
        labs( y = '1-FPR',x='') +
        ylim(0,1)+
        # scale_y_continuous(expand = c(0,0),limits = c(0, length0*1.55)) +
        ggtitle(figure_number,figure_name)+
        theme_classic()+
        theme(legend.position="none",
              axis.text.x=element_text(angle=45, hjust=1,face = 'plain',size=10),
              axis.text.y=element_text(face = 'plain',size=10),
              text = element_text(size=15,face = 'plain',family ="sans",colour = 'black'))+
        theme(legend.position="none",
              plot.title = element_text(hjust = -0.1, vjust=1,size=18))
      # plot.margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "cm"))
      
    }else{
      p <- ggplot(data_plot, aes(label, weight = num, fill = label)) +
        # geom_hline(yintercept = seq(10, 50, 10), color = 'gray') +
        geom_bar(color = "white", width = .65, position = 'dodge') +
        # scale_fill_manual(palette ="Set3")+  #values = c('#90BFF9','#05BE78','#FF0000','#B07FC0','#C5944E')
        # geom_errorbar(aes(ymin = num - se, ymax = num + se), width = 0.25, size = 0.3, position = position_dodge(0.7)) +
        labs( y = y_label0,x='') +
        ylim(0,1)+
        # scale_y_continuous(expand = c(0,0),limits = c(0, length0*1.55)) +
        ggtitle(figure_number,figure_name)+
        theme_classic()+
        theme(legend.position="none",
              axis.text.x=element_text(angle=45, hjust=1,face = 'plain',size=10),
              axis.text.y=element_text(face = 'plain',size=10),
              text = element_text(size=15,face = 'plain',family ="sans",colour = 'black'))+
        theme(legend.position="none",
              plot.title = element_text(hjust = -0.1, vjust=1,size=18))
      # plot.margin = margin(t = 0, r = 2, b = 0, l = 0, unit = "cm"))
    }
  }
  print(p)
  return(p)
}

# marker_split_by_cutoff ----
marker_split_by_cutoff <- function(Num0=5,
                                   cutoff_label='Train',
                                   model_method='lasso',
                                   class0='adjust'){
  if(model_method=='lasso'){
    load('./01_data/confonder_frame_all_lasso.RData')
  }else{
    load('./01_data/confonder_frame_all_rf.RData')
  }
  load('./01_data/feat_meta/All.feat.meta.nor.filter_taxon_sample.6.5.species.RData')
  feat_list0 <- my_data0$feat_list
  meta_list0 <- my_data0$meta_list
  feat_list0[c(1,2)] <- NULL
  meta_list0[c(1,2)] <- NULL
  load('./01_data/feat_meta/curatedMetagenomicData.nor.RData')
  feat_list0 <- c(my_data0$feat_list,feat_list0)
  meta_list0 <- c(my_data0$meta_list,meta_list0)
  feat_list0 <- my_pair.table(feat_list0)
  feat_list0$H2B_self <- NULL
  meta_list0$H2B_self <- NULL
  meta_list0 <- lapply(meta_list0, function(x){x[,c('Group','host_age','sex','BMI','disease_stage','country')]})
  
  # feat_list0$T2D_QinJ_2012 <- NULL
  # meta_list0$T2D_QinJ_2012 <- NULL
  
  # 1)disease PPI VS. non_PPI----
  feat_list <- feat_list0
  meta_list <- meta_list0
  for (i in names(feat_list)[-1]) {
    confonder_frame_all_t <- subset(confonder_frame_all,label==i)
    confonder_frame_all_t <- subset(confonder_frame_all_t,(Group==1))
    confonder_frame_all_t <- subset(confonder_frame_all_t,(Label==cutoff_label))
    # table(confonder_frame_all_t$pred_label)
    a <- intersect(confonder_frame_all_t$sample, rownames( feat_list[[i]]))
    if( length(a)==0){
      print('no way')
    }
    feat_list[[i]] <- feat_list[[i]][confonder_frame_all_t$sample,]
    meta_list[[i]] <- meta_list[[i]][confonder_frame_all_t$sample,]
    meta_list[[i]]$Group <- ifelse(confonder_frame_all_t$pred_label==1,'Case','Control')
    
    if(sum(table(meta_list[[i]]$Group)>=Num0)<2){
      feat_list[[i]] <- NULL
      meta_list[[i]] <- NULL
    }else{
      # print(rowSums(feat_list[[i]]==0))
      print(table(meta_list[[i]]$Group))
    }
    # print(table(meta_list[[i]]$Group))
  }
  marker <- my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,
                          is_plot=F,lda_cutoff=2,nproj_cutoff=1,level='species',change_name=F)
  
  save(marker,file = paste0('./01_data/model/new/Case_WGS_marker_PPI_VS_non-PPI','_',cutoff_label,
                            '_',model_method,'.RData'))
  # load('./01_data/model/new/Case_WGS_marker_PPI_VS_non-PPI_Test.RData')
  marker0 <- marker$marker_data
  marker1 <- subset(marker0,(nrproj>=2) & (class==class0))
  p <- marker_plot_heatmap(marker1,highlight=T)
  # p
  pdf(paste0('./02_figure/Case_PPI_WGS_marker_PPI_VS_non-PPI','_',cutoff_label,'_',model_method,'.pdf'),height = 30,width = 8)
  p
  dev.off()
  p1 <- p
  
  
  
  # 2)health VS. non_PPI disease----
  feat_list <- feat_list0
  meta_list <- meta_list0
  for (i in names(feat_list)[-1]) {
    confonder_frame_all_t <- subset(confonder_frame_all,label==i)
    confonder_frame_all_t <- subset(confonder_frame_all_t,(Group==-1)|((Group==1)&(pred_label==-1)))
    confonder_frame_all_t <- subset(confonder_frame_all_t,(Label==cutoff_label)) #Mean
    # table(confonder_frame_all_t$pred_label)
    a <- intersect(confonder_frame_all_t$sample, rownames( feat_list[[i]]))
    if( length(a)==0){
      print('no way')
    }
    feat_list[[i]] <- feat_list[[i]][confonder_frame_all_t$sample,]
    meta_list[[i]] <- meta_list[[i]][confonder_frame_all_t$sample,]
    meta_list[[i]]$Group <- ifelse(confonder_frame_all_t$Group==1,'Case','Control')
    # print(table(meta_list[[i]]$Group))
    if(sum(table(meta_list[[i]]$Group)>=Num0)<2){
      feat_list[[i]] <- NULL
      meta_list[[i]] <- NULL
    }else{
      # print(rowSums(feat_list[[i]]==0))
      print(table(meta_list[[i]]$Group))
    }
    # print(table(meta_list[[i]]$Group))
  }
  marker <- my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,
                          is_plot=F,lda_cutoff=2,nproj_cutoff=1,level='species',change_name=F)
  
  save(marker,file = paste0('./01_data/model/new/Case_WGS_marker_health_VS_non-PPI','_',cutoff_label,
                            '_',model_method,'.RData'))
  # load('./01_data/model/new/Case_WGS_marker_PPI_VS_non-PPI_Test.RData')
  marker0 <- marker$marker_data
  marker1 <- subset(marker0,(nrproj>=4) & (class==class0))
  p <- marker_plot_heatmap(marker1,highlight=T)
  # p
  pdf(paste0('./02_figure/Case_PPI_WGS_marker_health_VS_non-PPI','_',cutoff_label,'_',model_method,'.pdf'),height = 30,width = 8)
  p
  dev.off()
  p2 <- p
  
  # 2)health VS. PPI disease----
  feat_list <- feat_list0
  meta_list <- meta_list0
  for (i in names(feat_list)[-1]) {
    confonder_frame_all_t <- subset(confonder_frame_all,label==i)
    confonder_frame_all_t <- subset(confonder_frame_all_t,(Group==-1)|((Group==1)&(pred_label==1)))
    confonder_frame_all_t <- subset(confonder_frame_all_t,(Label==cutoff_label)) #Mean
    # table(confonder_frame_all_t$pred_label)
    a <- intersect(confonder_frame_all_t$sample, rownames( feat_list[[i]]))
    if( length(a)==0){
      print('no way')
    }
    feat_list[[i]] <- feat_list[[i]][confonder_frame_all_t$sample,]
    meta_list[[i]] <- meta_list[[i]][confonder_frame_all_t$sample,]
    meta_list[[i]]$Group <- ifelse(confonder_frame_all_t$Group==1,'Case','Control')
    # print(table(meta_list[[i]]$Group))
    if(sum(table(meta_list[[i]]$Group)>=Num0)<2){
      feat_list[[i]] <- NULL
      meta_list[[i]] <- NULL
    }else{
      # print(rowSums(feat_list[[i]]==0))
      print(table(meta_list[[i]]$Group))
    }
    # print(table(meta_list[[i]]$Group))
  }
  marker <- my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,
                          is_plot=F,lda_cutoff=2,nproj_cutoff=1,level='species',change_name=F)
  
  save(marker,file = paste0('./01_data/model/new/Case_WGS_marker_health_VS_PPI','_',cutoff_label,
                            '_',model_method,'.RData'))
  # load('./01_data/model/new/Case_WGS_marker_PPI_VS_non-PPI_Test.RData')
  marker0 <- marker$marker_data
  marker1 <- subset(marker0,(nrproj>=2) & (class==class0))
  p <- marker_plot_heatmap(marker1,highlight=T)
  # p
  pdf(paste0('./02_figure/Case_PPI_WGS_marker_health_VS_PPI','_',cutoff_label,'_',model_method,'.pdf'),height = 30,width = 8)
  p
  dev.off()
  p3 <- p
  
  # health_enrich <- marker0$scientific_name[(as.numeric(marker0$LDA)<0)]
  # disease_enrich <- marker0$scientific_name[(as.numeric(marker0$LDA)>0)]
  # marker0 <- subset(marker0,scientific_name %in% health_enrich)
  # marker0 <- subset(marker0,!scientific_name %in% disease_enrich)
  P <- list(p1=p1,p2=p2,p3=p3)
  return(P)
}

# model_plot_ROC ----
# col is sample
#' @param row_t plot which row of feat_meta_info; only when is.null(model.path)=T can work
#' @param model_method model method
#' @param self_project which cohort is selected to be examined
#' @param cutoff AUC cutoff for plot
#' @param model.path model path


model_plot_ROC <- function(row_t=1,model_method='lasso',self_project='PPI_self',cutoff=0,model.path=NULL){
  # for (row_t in seq_len(nrow(feat_meta_info))){
  # library(tidyr)
  # library(dplyr)
  if(is.null(model.path)){
    Group_t <- feat_meta_info$Group[row_t]
    PPI_label_t <- feat_meta_info$PPI_label[row_t]
    data_type_t <- feat_meta_info$data_type[row_t]
    feat_meta_path_t <- feat_meta_info$feat_meta_path[row_t]
    print(paste(Group_t,data_type_t,PPI_label_t,sep = '_'))
    print(feat_meta_path_t)
    is.plot <- paste0('02_figure/',model_method,'/ROC_plot_',paste(Group_t,data_type_t,PPI_label_t,sep = '_'),'.pdf')
    model.adj <- get(load(paste0('D:/lm_project/PPI_cross/01_data/largeRData/model/',model_method,'/model.',
                                 paste(row_t,Group_t,data_type_t,PPI_label_t,sep = '_'),'.only_external.RData')))
  }else{
    model.adj <- get(load(model.path))
    is.plot <- paste0(model.path,'.ROC.pdf')
  }
  
  
  # ROC
  ROCdata <- data.frame()
  for (pro_t in rownames((model.adj[["result"]]))) {
    if(pro_t==self_project){
      AUC_t <- model.adj$result[self_project,pro_t]
      model_siamcat <- model.adj$models_top[[self_project]]
      models_Y <- model_siamcat$model_back
      sensitivities <- models_Y@eval_data$roc$sensitivities
      specificities <- models_Y@eval_data$roc$specificities
      ROCdata_t <- data.frame(x=1-specificities,y=sensitivities,Group=paste0(pro_t,' AUC: ',AUC_t),AUC=AUC_t)
      ROCdata_t <- arrange(ROCdata_t,y,x)
      ROCdata <- rbind(ROCdata,ROCdata_t)
    }else{
      AUC_t <- model.adj$result[self_project,pro_t]
      model_siamcat <- model.adj$cross_models[[paste0(self_project,'-',pro_t)]]
      models_Y <- model_siamcat$model_back
      sensitivities <- models_Y@eval_data$roc$sensitivities
      specificities <- models_Y@eval_data$roc$specificities
      ROCdata_t <- data.frame(x=1-specificities,y=sensitivities,Group=paste0(pro_t,' AUC: ',AUC_t),AUC=AUC_t)
      ROCdata_t <- arrange(ROCdata_t,y,x)
      ROCdata <- rbind(ROCdata,ROCdata_t)
    }
    
  }
  
  ROCdata <- subset(ROCdata,AUC>=cutoff)
  
  
  # library(ggplot2)
  # library(ROCR)
  # library(ggsci)
  # library(dplyr)
  # if(!is.null(is.plot)){
  p1 <- ROCdata%>%
    ggplot(aes(x=x, y=y)) + 
    # geom_line(aes(color=cut(y, c(-1,0.5,0.7,0.9,1)))) +
    # geom_point(aes(color=cut(y, c(-1,0.5,0.7,0.9,1)))) +
    geom_line(aes(color=Group)) +
    # geom_point() +
    theme_bw()+
    labs(color="Group", x="False positive rate", y="True positive rate") +
    # scale_color_brewer(palette = "Reds") +
    scale_color_igv()+
    theme(legend.position = c(0.7,0.3),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 10))
  # geom_text(aes(x=1, y=0, label=paste0("Mean AUC: ",'mean_AUC','; Accuracy: ','AUT_type2')),
  # hjust="right", vjust="bottom",size=5)
  print(p1)
  
  ggsave(filename = is.plot,width = 5,height = 5)
  # }
  # }
  return(ROCdata)
}

# model_plot_heatmap ----
# col is sample
#' @param row_t plot which row of feat_meta_info; only when is.null(model.path)=T can work
#' @param self_project which cohort is selected to be examined
#' @param model.path model path

model_plot_heatmap <- function(row_t=1,model_method='lasso',self_project='PPI_self',model.path=NULL){
  # for (row_t in seq_len(nrow(feat_meta_info))){
  
  if(is.null(model.path)){
    Group_t <- feat_meta_info$Group[row_t]
    PPI_label_t <- feat_meta_info$PPI_label[row_t]
    data_type_t <- feat_meta_info$data_type[row_t]
    feat_meta_path_t <- feat_meta_info$feat_meta_path[row_t]
    print(paste(Group_t,data_type_t,PPI_label_t,sep = '_'))
    print(feat_meta_path_t)
    is.plot <- paste0('02_figure/',model_method,'/Heatmap_plot_',paste(Group_t,data_type_t,PPI_label_t,sep = '_'),'.pdf')
    model.adj <- get(load(paste0('D:/lm_project/PPI_cross/01_data/largeRData/model/',model_method,'/model.',
                                 paste(row_t,Group_t,data_type_t,PPI_label_t,sep = '_'),'.RData')))
  }else{
    model.adj <- get(load(model.path))
    is.plot <- paste0(model.path,'.hp.pdf')
  }
  
  # heatmap
  hp <- my_auc.heatmap(datas=model.adj$result,models=paste0(model_method,'_adj_nc'),excluded=F,excluded_name=c('PPI_self'))
  pdf(is.plot,
      width = length(model.adj$models_top)+2,
      height = length(model.adj$models_top)-1)
  hp
  dev.off()
  
  return(hp)
}


# model_plot_ROC_and_heatmap ----
# col is sample
#' @param row_t plot which row of feat_meta_info; only when is.null(model.path)=T can work
#' @param self_project which cohort is selected to be examined
#' @param model.path model path


model_plot_ROC_and_heatmap <- function(row_t=1,model_method='lasso',self_project='PPI_self',model.path=NULL){
  # for (row_t in seq_len(nrow(feat_meta_info))){
  
  if(is.null(model.path)){
    Group_t <- feat_meta_info$Group[row_t]
    PPI_label_t <- feat_meta_info$PPI_label[row_t]
    data_type_t <- feat_meta_info$data_type[row_t]
    feat_meta_path_t <- feat_meta_info$feat_meta_path[row_t]
    print(paste(Group_t,data_type_t,PPI_label_t,sep = '_'))
    print(feat_meta_path_t)
    is.plot <- paste0('02_figure/',model_method,'/ROC_plot_',paste(Group_t,data_type_t,PPI_label_t,sep = '_'),'.pdf')
    is.plot.hp <- paste0('02_figure/',model_method,'/Heatmap_plot_',paste(Group_t,data_type_t,PPI_label_t,sep = '_'),'.pdf')
    model.adj <- get(load(paste0('D:/lm_project/PPI_cross/01_data/largeRData/model/',model_method,'/model.',
                                 paste(row_t,Group_t,data_type_t,PPI_label_t,sep = '_'),'.only_external.RData')))
  }else{
    model.adj <- get(load(model.path))
    is.plot <- paste0(model.path,'.ROC.pdf')
    is.plot.hp <- paste0(model.path,'.hp.pdf')
  }
  
  # heatmap
  hp <- my_auc.heatmap(datas=model.adj$result,models=paste0(model_method,'_adj_nc'),excluded=F,excluded_name=c('PPI_self'))
  pdf(is.plot.hp,
      width = length(model.adj$models_top)+2,
      height = length(model.adj$models_top)-1)
  hp
  dev.off()
  
  # ROC
  ROCdata <- data.frame()
  for (pro_t in names(model.adj$models_top)) {
    if(pro_t==self_project){
      AUC_t <- model.adj$result[self_project,pro_t]
      model_siamcat <- model.adj$models_top[[self_project]]
      models_Y <- model_siamcat$model_back
      sensitivities <- models_Y@eval_data$roc$sensitivities
      specificities <- models_Y@eval_data$roc$specificities
      ROCdata_t <- data.frame(x=1-specificities,y=sensitivities,Group=paste0(pro_t,' AUC: ',AUC_t))
      ROCdata_t <- arrange(ROCdata_t,y,x)
      ROCdata <- rbind(ROCdata,ROCdata_t)
    }else{
      AUC_t <- model.adj$result[self_project,pro_t]
      model_siamcat <- model.adj$cross_models[[paste0(self_project,'-',pro_t)]]
      models_Y <- model_siamcat$model_back
      sensitivities <- models_Y@eval_data$roc$sensitivities
      specificities <- models_Y@eval_data$roc$specificities
      ROCdata_t <- data.frame(x=1-specificities,y=sensitivities,Group=paste0(pro_t,' AUC: ',AUC_t))
      ROCdata_t <- arrange(ROCdata_t,y,x)
      ROCdata <- rbind(ROCdata,ROCdata_t)
    }
    
  }
  # library(ggplot2)
  # library(ROCR)
  # library(ggsci)
  # library(dplyr)
  # if(!is.null(is.plot)){
  p1 <- ROCdata%>%
    ggplot(aes(x=x, y=y)) + 
    # geom_line(aes(color=cut(y, c(-1,0.5,0.7,0.9,1)))) +
    # geom_point(aes(color=cut(y, c(-1,0.5,0.7,0.9,1)))) +
    geom_line(aes(color=Group)) +
    # geom_point() +
    theme_bw()+
    labs(color="Group", x="False positive rate", y="True positive rate") +
    # scale_color_brewer(palette = "Reds") +
    scale_color_igv()+
    theme(legend.position = c(0.7,0.3),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 10))
  # geom_text(aes(x=1, y=0, label=paste0("Mean AUC: ",'mean_AUC','; Accuracy: ','AUT_type2')),
  # hjust="right", vjust="bottom",size=5)
  print(p1)
  
  ggsave(filename = is.plot,width = 5,height = 5)
  # }
  # }
  return()
}
# heatmap correlation ----
# row is feature
# col is sample

heatmap.correlation <- function(F1,F2,p_adjust=T){
  
  ## Spearman
  x4<-data.frame(cor=NULL)
  x5<-data.frame(cor=NULL)
  for (i0 in rownames(F1)){
    for (j0  in rownames(F2)){
      x4[i0,j0] <- as.numeric(cor(t(F1[i0,]),t(F2[j0,]),method = "spearman"))
      x5[i0,j0] <- as.numeric(cor.test(t(F1[i0,]),t(F2[j0,]),method = "spearman")[["p.value"]])
    }
  }
  x5_non_adjust <- x5
  
  
  x5_adjust <- c()
  for (row_t in seq_len(nrow(x5))) {
    x5_adjust <- c(x5_adjust,x5[row_t,])
  }
  x5_adjust <- p.adjust(x5_adjust,method = "BH")
  for (row_t in seq_len(nrow(x5))) {
    c_from <- (row_t*ncol(x5)-(ncol(x5)-1))
    c_end <- row_t*ncol(x5)
    x5[row_t,] <- x5_adjust[c_from:c_end]
  }
  
  
  
  
  x6 <- x5
  x6[x6>0.05] <- 0
  x6 <- x6[,colSums(x6)!=0]
  x4 <- x4[,colnames(x6)]
  x5 <- x5[,colnames(x6)]

  
  
  data_mark <- data.frame()
  for(i in 1:nrow(x5)){
    for(j in 1:ncol(x5)){
      if(x5[i,j] <= 0.001)
      {
        data_mark[i,j]="***"
      }
      else if(x5[i,j] <= 0.01 && x5[i,j] > 0.001)
      {
        data_mark[i,j]="**"
      }
      else if(x5[i,j] <= 0.05 && x5[i,j] > 0.01)
      {
        data_mark[i,j]="*"
      }
      else
      {
        data_mark[i,j]=""
      }
    }
  }
  
  data_mark_non_adjust <- data.frame()
  for(i in 1:nrow(x5_non_adjust)){
    for(j in 1:ncol(x5_non_adjust)){
      if(x5_non_adjust[i,j] <= 0.001)
      {
        data_mark_non_adjust[i,j]="***"
      }
      else if(x5_non_adjust[i,j] <= 0.01 && x5_non_adjust[i,j] > 0.001)
      {
        data_mark_non_adjust[i,j]="**"
      }
      else if(x5_non_adjust[i,j] <= 0.05 && x5_non_adjust[i,j] > 0.01)
      {
        data_mark_non_adjust[i,j]="*"
      }
      else
      {
        data_mark_non_adjust[i,j]=""
      }
    }
  }
  # annotation_row <- read.table('function_index_annotation.csv',sep = ',',header = T,row.names = 1)
  # library(pheatmap)
  # library(ggplotify)
  #navy 
  if(p_adjust==T){
    p <- pheatmap(x4, color = colorRampPalette(c("#7BA7A6", "white", "firebrick3"))(50),display_numbers = data_mark,fontsize_number=16,cluster_row = F,angle_col = "45",cellwidth = 29, cellheight = 13,silent=T)
  }else{
    p <- pheatmap(x4, color = colorRampPalette(c("#7BA7A6", "white", "firebrick3"))(50),display_numbers = data_mark_non_adjust,fontsize_number=16,cluster_row = F,angle_col = "45",cellwidth = 29, cellheight = 13,silent=T)
  }
  
  p <- as.ggplot(p)
  correlation_result <- list(R=x4,p_value=x5,p=p)
  return(correlation_result)
}

# PCA ----
PCA <- function(meta,
                feat,
                title='PCA bray distance',
                method0='bray',
                x0=0.3,
                y0=-0.2,
                only.p1=F
){
  # feat col is feature, row is sample
  # meta row is sample
  # library('vegan')
  meta <- meta[rownames(feat),]
  # Matriax
  data.dist=vegdist(feat[,c(1,2)],method0)
  
  # message(paste0('You can use ',paste(c( "euclidean", "bray", "jaccard", "manhattan", "canberra", "clark",  "kulczynski", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq" ,"chord"),collapse = ',')))
  
  # ADONIS/Permannova
  # adonis.result <- adonis(feat ~ Group,
  #                         data = meta,
  #                         distance = "bray",
  #                         permutations = 999)
  # adonis.result$aov.tab$Pr[1]
  
  # Anosim 
  anosim.result<-anosim(data.dist,meta$Group,permutations = 999)
  # summary(anosim.result)
  
  # PCA
  # library(ape)
  pcoa_ind = pcoa(data.dist)
  pc = as.data.frame(pcoa_ind$vectors)
  colnames(pc)=gsub('Axis.','PCo',colnames(pc))
  pc$Group=as.factor(meta$Group)
  color=c( "firebrick3","navy","#91D1C2B2","#3C5488B2" ,
           "#8491B4B2", "#DC0000B2", 
           "#7E6148B2","yellow", 
           "darkolivegreen1", "lightskyblue", 
           "darkgreen", "deeppink", "khaki2", 
           "firebrick", "brown1", "darkorange1", 
           "cyan1", "royalblue4", "darksalmon", 
           "darkgoldenrod1", "darkseagreen", "darkorchid")
  
  # Plot
  p1 <- ggplot(pc,aes(PCo1,PCo2,color=Group,Group = Group))+
    geom_hline(yintercept=0,linetype=2,color='black')+
    geom_vline(xintercept=0,linetype=2,color='black')+ 
    stat_ellipse(size=0.5,level = 0.95)+ # level meaning confidence
    # ggtitle(paste0(project[i,c('project_id','disease',"data_type",'taxon')],collapse = '_'))+ 
    geom_point(aes(shape=Group),size=2.5)+ 
    annotate('text',x=x0,y=y0,label = paste0('R = ',round(anosim.result[["statistic"]],2), ',','p = ',round(anosim.result[["signif"]],4)))+
    scale_color_manual(values =color[1:length(unique(meta$Group))]) +
    scale_fill_manual(values = color[1:length(unique(meta$Group))]) +
    theme_bw()+
    theme(panel.grid=element_blank(),axis.text = element_text(size=12),axis.title = element_text(size=14)) 
  
  if(!is.null(title)){
    p1 <- p1+ggtitle(title)
  }
  
  if(only.p1==T){
    p=p1
  }else{
    
    # library('dplyr')
    # library('ggpubr')
    
    # 写一个生成posi.x，和 posi.y的函数
    # posi.x=0.2
    # posi.y=0.5
    
    merged <- pc[,c('PCo1','PCo2','Group')]
    
    length.out0 <- length(unique(meta$Group))
    length.out0 <- length.out0*(length.out0-1)/2
    stat.test <- compare_means(
      PCo2~Group,data = merged,
      method = "t.test"
    )%>%mutate(y.position = seq(max(merged$PCo2),max(merged$PCo2)+max(merged$PCo2)/2,length.out=length.out0))
    x=stat.test$p.adj
    stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
    p2 <- ggboxplot(merged,x= 'Group', y="PCo2",fill="Group",palette = color[1:length(unique(meta$Group))])+
      # clean_theme()+
      xlab('')+ylab('')+
      # ylim(-0.27,0.27)+
      theme(legend.position = 'none',axis.text.x = element_blank())
    p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
    # p2
    
    length.out0 <- length(unique(meta$Group))
    length.out0 <- length.out0*(length.out0-1)/2
    stat.test <- compare_means(
      PCo1~Group,data = merged,
      method = "t.test"
    )%>%mutate(y.position = seq(max(merged$PCo1),max(merged$PCo1)+max(merged$PCo1)/2,length.out=length.out0))
    x=stat.test$p.adj
    stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
    p3 <- ggboxplot(merged,x= 'Group', y="PCo1",fill="Group",palette = color[1:length(unique(meta$Group))])+
      # clean_theme()+
      xlab('')+ylab('')+
      # ylim(-0.3,0.75)+
      theme(legend.position = 'none',axis.text.y = element_blank())+
      ggpubr::rotate()
    p3 <- p3+stat_pvalue_manual(stat.test,label = "p.adj.signif")
    # p3
    
    p=ggarrange(p1, p2, p3, NULL, ncol = 2, nrow = 2, align = "hv", widths = c(3, 1),
                heights = c(3,2),common.legend = T)
  }
  
  # ggsave('905.pdf',height=12,width = 18,units = cm)
  return(p)
}


# kmeans_pca_plot ----
kmeans_pca_plot <- function(data, centers,centers0) {
  # library(ggplot2)
  # library(cluster)
  # 储存不同聚类个数的总内部离差平方和（WCSS）
  wcss <- vector("numeric", length = centers)  # 设置聚类个数的范围
  
  # 计算不同聚类个数的总内部离差平方和
  for (k in 1:centers) {
    kmeans_result <- kmeans(data, centers = k)
    wcss[k] <- kmeans_result$tot.withinss
  }
  
  # 绘制总内部离差平方和与聚类个数的折线图
  plot_data <- data.frame(K = 1:centers, WCSS = wcss)
  elbow_plot <- ggplot(plot_data, aes(x = K, y = WCSS)) +
    geom_line() +
    geom_point() +
    labs(title = "Elbow Method", x = "Number of Clusters (K)", y = "WCSS")
  
  # 执行 K 均值聚类
  kmeans_result <- kmeans(data, centers = centers0)
  
  # 将聚类结果添加到原始数据集
  clustered_data <- cbind(data, cluster = as.factor(kmeans_result$cluster),Group=as.factor(meta$Group))
  
  # 创建 PCA 对象
  pca_result <- prcomp(data)
  
  # 提取前两个主成分的得分
  pca_scores <- data.frame(pca_result$x[, c(1, 2)])
  
  # 组合 PCA 得分和聚类结果
  plot_data <- cbind(pca_scores, clustered_data)
  
  # 绘制 PCA 坐标图，并按聚类结果着色
  pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster,shape = Group)) +
    geom_point() +
    labs(title = "K-means Clustering with PCA", x = "PC1", y = "PC2")+
    theme_bw()
  
  return(list(elbow_plot=elbow_plot, pca_plot=pca_plot))
}

# network_plot -----

network_plot <- function(project_id='PPI_self',r_cutoff=0.6,isplot=T){
  
  # 导入所需的包
  # library(igraph)
  # library(psych)
  
  # 读取菌群丰度表
  my_data0 <- get(load(paste0('01_data/feat_meta/standard/data_disease_WGS_PPI_unknown.RData')))
  feat_list <- my_data0$feat_list
  meta_list <- my_data0$meta_list
  feat <- feat_list[[project_id]]
  meta <- meta_list[[project_id]]
  feat <- feat[,colSums(feat>0)>2]
  abundance <- feat
  if(project_id=='PPI_self'){
    load('01_data/PPI_species_df.RData')
    PPI_species <- subset(PPI_species_df,PPI_species_df$type=='Marker species')$node
  }else{
    load('01_data/PPI_species.RData')
  }
  
  
  
  # abundance <- abundance[,PPI_species]
  
  
  
  
  type <- get(load('01_data/PPI_species_df.RData'))
  cor <- corr.test(abundance, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=FALSE) # ci=FALSE,不进行置信区间计算，数据量较大时，可以加快计算速度。
  cor.r <- data.frame(cor$r) # 提取R值
  cor.p <- data.frame(cor$p) # 提取p值
  all_sp <- colnames(abundance)
  
  
  if(length(all_sp[all_sp%in%type$node])!=0){
    type <- rbind(subset(type, node %in% all_sp),
                  data.frame('node'=all_sp[!all_sp%in%type$node],
                             'type'='Species'))
  }
  
  
  
  colnames(cor.r) = rownames(cor.r)
  colnames(cor.p) = rownames(cor.p) # 变量名称中存在特殊字符，为了防止矩阵行名与列名不一致，必须运行此代码。
  # write.csv(cor.r,paste0("03_table/network/",project_id,".cor.r.csv"),quote = FALSE,row.names = TRUE) # 保存结果到本地
  # write.csv(cor.p,paste0("03_table/network/",project_id,".cor.p.csv"),quote = FALSE,row.names = TRUE)
  

  ### 2.2 确定相关性关系 保留p<=0.05且abs(r)>=0.6的变量间相关关系
  # cor.r[abs(cor.r) < 0.6 | cor.p > 0.05] = 0
  # cor.r[cor.p > 0.05] = 0
  # cor.r = as.matrix(cor.r)
  # g = graph_from_adjacency_matrix(as.matrix(cor.r),mode = "undirected",weighted = TRUE,diag = FALSE)
  
  #### 将数据转换为long format进行过滤，后面绘制网络图需要节点和链接数据，在这一步可以完成格式整理
  cor.r$node1 = rownames(cor.r) 
  cor.p$node1 = rownames(cor.p)
  # library(tidyr)
  r = cor.r %>% 
    gather(key = "node2", value = "r", -node1) %>%
    data.frame()
  # r <- subset(r,(node1%in%PPI_species)|(node2%in%PPI_species))
  r <- subset(r,(node1%in%PPI_species))
  p = cor.p %>% 
    gather(key = "node2", value = "p", -node1) %>%
    data.frame()
  # p <- subset(p,(node1%in%PPI_species)|(node2%in%PPI_species))
  p <- subset(p,(node1%in%PPI_species))
  # head(r)
  # head(p)
  
  #### 将r和p值合并为一个数据表
  cor.data <- merge(r,p,by=c("node1","node2"))
  # cor.data
  # library(dplyr)
  #### 保留p<=0.05且abs(r)>=0.6的变量间相关关系，并添加网络属性
  cor.data <- cor.data %>%
    filter(abs(r) >= r_cutoff, p <= 0.05, node1 != node2) %>%  #abs(r) >= 0.6,
    mutate(
      linetype = ifelse(r > 0,"positive","negative"), # 设置链接线属性，可用于设置线型和颜色。
      linesize = abs(r) # 设置链接线宽度。
    ) # 此输出仍有重复链接，后面需进一步去除。
  # head(cor.data)
  
  
  if(!isTRUE(isplot)){
    return(cor.data)
  }else{
    ### 3.1 网络图节点属性整理
    #### 计算每个节点具有的链接数
    c(as.character(cor.data$node1),as.character(cor.data$node2)) %>%
      as_tibble() %>%
      group_by(value) %>%
      summarize(n=n()) -> vertices
    colnames(vertices) <- c("node", "n")
    # head(vertices)
    
    #### 添加变量分类属性
    vertices <- vertices %>%
      select(-n) %>% # 因为此处的n不准确，所以先删除。
      left_join(type,by="node")
    
    #### 网络图中节点会按照节点属性文件的顺序依次绘制，为了使同类型变量位置靠近，按照节点属性对节点进行排序。
    vertices$type = factor(vertices$type,levels = unique(vertices$type))
    vertices = arrange(vertices,type)
    # write.csv(vertices,paste0("03_table/network/",project_id,".vertices.csv"),quote = FALSE,
    #           row.names = FALSE)
    # head(vertices)
    
    
    ### 3.2 构建graph结构数据
    g <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
    # g
    # vcount(g) # 节点数目：82
    # ecount(g) # 链接数:297
    #get.vertex.attribute(g) # 查看graph中包含的节点属性
    #get.edge.attribute(g) # 查看graph中包含的链接属性
    
    ### 3.3 简单图
    # is.simple(g) # 非简单图，链接数会偏高，所以需要转换为简单图。
    E(g)$weight <- 1
    g <- igraph::simplify(g,
                          remove.multiple = TRUE,
                          remove.loops = TRUE,
                          edge.attr.comb = "first")
    g <- delete.vertices(g,which(degree(g) == 0)) # 删除孤立点
    # is.simple(g)
    E(g)$weight <- 1
    # is.weighted(g)
    # vcount(g) # 节点数目：21
    # ecount(g) # 链接数:33
    
    
    ### 3.4 计算节点链接数
    V(g)$degree <- degree(g)
    #vertex.attributes(g)
    #edge.attributes(g) 
    # g
    
    ### 3.5 graph保存到本地
    # write.graph(g,file = paste0("02_figure/network/",project_id, ".all.gml"),format="gml") # 直接保存graph结构，gml能保存的graph信息最多。
    net.data  <- igraph::as_data_frame(g, what = "both")$edges # 提取链接属性
    # write.csv(net.data,file=paste0("03_table/network/",project_id,".net.data.csv"),quote = FALSE,row.names = FALSE) # 保存结果到本地。
    # head(net.data)
    
    vertices  <- igraph::as_data_frame(g, what = "both")$vertices # 提取节点属性
    # write.csv(vertices,paste0("03_table/network/",project_id,".vertices.csv"),quote = FALSE,
    #           row.names = FALSE)
    # head(vertices) # 直接读入之前保存的链接和节点属性文件，之后可直接生成graph或用于其他绘图软件绘图。
    
    
    ### 4.1 准备网络图布局数据
    #?layout_in_circle # 帮助信息中，有其它布局函数。
    layout1 <- layout_in_circle(g) # 径向布局适合节点较少的数据。
    layout2 <- layout_with_fr(g) # fr布局。
    layout3 <- layout_on_grid(g) # grid布局。
    # head(layout1)
    
    ### 4.2 设置绘图颜色
    
    color <- c("#41B3C2",
               "#C9D8C5",
               "#FFE7AE")
    # "#FF9B00"
    
    # library(ggsci)
    # color <- pal_aaas(alpha = 1)(3) #这里是指提取前10个颜色
    # library(scales)
    # show_col(color)
    
    names(color) <- unique(V(g)$type) # 将颜色以节点分类属性命名
    # names(color) <- c("Correlation Species","Species","Marker species") # 将颜色以节点分类属性命名
    V(g)$point.col <- color[match(V(g)$type,names(color))] # 设置节点颜色。
    #names(color2) <- unique(V(g)$type) # 如果想要节点颜色与背景颜色不一致，则可以为节点单独设置一个颜色集。
    #V(g)$point.col <- color2[match(V(g)$type,names(color2))]
    
    #### 边颜色按照相关性正负设置
    #E(g)$color <- ifelse(E(g)$linetype == "positive",rgb(255,215,0,maxColorValue = 255),"gray50")
    E(g)$color <- ifelse(E(g)$linetype == "positive","red",rgb(0,147,0,maxColorValue = 255))
    
    
    
    ### 4.3 绘制径向布局网络图
    pdf(paste0("02_figure/network/",project_id,".network_group_circle.pdf"),family = "Times",width = 10,height = 12)
    par(mar=c(5,2,1,2))
    plot.igraph(g, layout=layout1,#更多参数设置信息?plot.igraph查看。
                
                ##节点颜色设置参数##
                vertex.color=V(g)$point.col,
                vertex.frame.color ="black",
                vertex.border=V(g)$point.col,
                ##节点大小设置参数##
                vertex.size=V(g)$degree*1.8,
                ##节点标签设置参数##
                vertex.label=g$name,
                vertex.label.cex=0.8,
                #vertex.label.dist=0, # 标签距离节点中心的距离，0表示标签在节点中心。
                #vertex.label.degree = 0, # 标签对于节点的位置，0-右，pi-左，-pi/2-上，pi/2-下。
                vertex.label.col="black",
                
                ##链接属性参数以edge*开头##
                edge.arrow.size=0.5,
                edge.width=abs(E(g)$r)^2*5,
                edge.curved = TRUE
    )
    ##设置图例，与plot.igraph()函数一起运行##
    legend(
      title = "Type",
      list(x = min(layout1[,1])-0.2,
           y = min(layout1[,2])-0.17), # 图例的位置需要根据自己的数据进行调整，后面需要使用AI手动调整。
      legend = c(unique(V(g)$type)),
      fill = color,
      #pch=1
    )
    
    legend(
      title = "|r-value|",
      list(x = min(layout1[,1])+0.4,
           y = min(layout1[,2])-0.17),
      legend = c(0.6,0.8,1.0),
      col = "black",
      lty=1,
      lwd=c(0.6,0.8,1.0)^2*5,
    )
    legend(
      title = "Correlation (±)",
      list(x = min(layout1[,1])+0.8,
           y = min(layout1[,2])-0.17),
      legend = c("positive","negative"),
      col = c("red",rgb(0,147,0,maxColorValue = 255)),
      lty=1,
      lwd=1
    )
    legend(
      title = "Degree",
      list(x = min(layout1[,1])+1.2,
           y = min(layout1[,2])-0.17),
      legend = c(1,seq(0,8,2)[-1]),# V(g)$degree。
      col = "black",
      pch=1,
      pt.lwd=1,
      pt.cex=c(1,seq(0,8,2)[-1]) 
      # 设置相同数值，igraph与legend生成的节点大小不一样，
      # 图例节点差不多是图节点的2倍大，这里通过设置图节点大小为degree*2解决，需要手动调整。
    )
    dev.off()
    return(cor.data)
  }
  
  
}


# marker_times_folds ----
marker_times_folds <- function(feat_list,meta_list,
                       check.con.before=NULL,check.con.after=NULL,
                       add_group=T,is_combat=T,re_scale=T,
                       is_plot=F,lda_cutoff=2,nproj_cutoff=1,
                       level='genus',change_name=F,do.con=T,times0=3,folds0=5){
  # library(microbiomeMarker)
  # library(phyloseq)
  # library(ggplot2)
  # library(RMySQL)
  # library(tidyr)
  # library(dplyr)
  split_list <- list()
  
    
  for (project_t in names(feat_list)) {
    split_list[[project_t]] <- list()
    
    
    feat_t <- feat_list[[project_t]]
    meta_t <- meta_list[[project_t]]
    case_rows <- rownames(meta_t[meta_t$Group == "Case", ])
    control_rows <- rownames(meta_t[meta_t$Group == "Control", ])
    
    for (times_t in 1:times0) {
      split_list[[project_t]][[times_t]] <- list()
      
      vec <- case_rows
      base_size <- floor(length(vec) / folds0)
      splits <- sample(rep(1:folds0, each = base_size, length.out = base_size*folds0))
      subsets_case <- split(vec, splits,drop = FALSE)
      
      vec <- control_rows
      base_size <- floor(length(vec) / folds0)
      splits <- sample(rep(1:folds0, each = base_size, length.out = base_size*folds0))
      subsets_control <- split(vec, splits,drop = FALSE)
      
      subsets <- Map(c,subsets_case,subsets_control)
      
      for (folds_t in 1:folds0) {
        split_list[[project_t]][[times_t]][[folds_t]] <- list()
        split_list[[project_t]][[times_t]][[folds_t]] <- unlist(subsets[-folds_t]) %>% as.character()
      }
    }
  }
  
  marker_list <- list()
  for (times_t in 1:times0) {
    marker_list[[times_t]] <- list()
    for (folds_t in 1:folds0) {
      marker_list[[times_t]][[folds_t]] <- list()
      
      feat_list_name <- lapply(split_list, function(x){
        x[[times_t]][[folds_t]] 
      })
      feat_list_t <- list()
      meta_list_t <- list()
      for (project_t in names(feat_list)) {
        feat_list_t[[project_t]] <- feat_list[[project_t]][feat_list_name[[project_t]],]
        meta_list_t[[project_t]] <- meta_list[[project_t]][feat_list_name[[project_t]],]
      }
      
      m0=my_marker.plot_new(markers_data=NULL,feat_list_t,meta_list_t,
                            lda_cutoff,nproj_cutoff,level,change_name,cut=F)
      if(do.con==T){
        feat_all_adj=list()
        if(!is.null(check.con.before)){
          check.con.before=paste0(check.con.before,names(feat_list_t),'.pdf')
        }
        if(!is.null(check.con.after)){
          check.con.after=paste0(check.con.after,names(feat_list_t),'.pdf')
        }
        
        for(i in names(feat_list_t)){
          feat_all_adj[[i]]=my_adj(feat_list_t[[i]],meta_list_t[[i]],check.con.before[i],check.con.after[i],add_group,is_combat,re_scale)
        }
        feat=lapply(feat_all_adj,function(x) x$feat_new)
        confounder=lapply(feat_all_adj,function(x) x$confounder)
        meta=lapply(feat_all_adj,function(x) x$meta_new)
        
        for (p_t_i in names(confounder)) {
          if(is.null(confounder[[p_t_i]])){
            print(p_t_i)
            feat[[p_t_i]] <- feat_list_t[[p_t_i]]
            meta[[p_t_i]] <- meta_list_t[[p_t_i]]
            
          }
        }
        m1=my_marker.plot_new(markers_data=NULL,feat_list=feat,meta_list=meta,
                              lda_cutoff,nproj_cutoff,level,change_name,cut=F)
        
        m0$class='origin'
        m1$class='adjust'
        m=rbind(m0,m1)
      }else{
        m1=m0
        m0$class='origin'
        m1$class='adjust'
        m=rbind(m0,m1)
      }
      
      if(is_plot==T){
        m1$project_id=paste0(m1$project_id,'_adj')
        mm=rbind(m0,m1)
        mm$project_id=factor(mm$project_id,levels = c(unique(m0$project_id),paste0(unique(m0$project_id),'_adj')))
        mm=mm[order(mm$project_id),]
        mm=my_marker.plot_new(mm,lda_cutoff = lda_cutoff,
                              nproj_cutoff = nproj_cutoff,level = level)
      }
      
      marker_list[[times_t]][[folds_t]] <- list(marker_data=m,confounder=confounder) 
    }
  }
  list(marker_list=marker_list,split_list=split_list)
}

# pre_prossess_function ----
pre_prossess_function <- function(genus.data,filter.taxon=T,filter.sample=T){
  genus.data <- t(genus.data)
  
  # #exclude sample prevelence < 2
  # genus.data <- genus.data[,colSums(genus.data!=0)>2]
  # #exclude taxon prevelence < 2
  # genus.data <- genus.data[rowSums(genus.data!=0)>2,]
  # #normalization
  genus.data <- genus.data/matrix(colSums(genus.data), nrow = nrow(genus.data) , ncol = ncol(genus.data) , byrow=TRUE)
  # # exclude abundence < 0.001
  # genus.data <- genus.data[rowSums(genus.data>0.001)>0,]
  #
  if(filter.taxon==T){
    ori_len <- nrow(genus.data)
    #exclude prevelence < 2
    genus.data <- genus.data[rowSums(genus.data!=0)>2,]
    #exclude all abundance < 0.001
    genus.data <- genus.data[rowSums(genus.data>=0.001)>0,]
    now_len <- nrow(genus.data)
    print(paste0('++++ filter ',ori_len-now_len,' taxon!'))
    print(paste0('++++ retain ',now_len,' taxon!'))
    
  }
  if(filter.sample==T){
    #exclude prevelence < 2
    ori_len <- ncol(genus.data)
    genus.data <- genus.data[,colSums(genus.data!=0)>2]
    now_len <- ncol(genus.data)
    print(paste0('++++ filter ',ori_len-now_len,' samples!'))
    print(paste0('++++ retain ',now_len,' samples!'))
  }
  
  
  genus.data <- t(genus.data)
  return(genus.data)
}






# metaphlan4.species.name ----
#row is species
metaphlan4.species.name <- function(data,name_exclude=T){
  data_species <- data[-grep('t__',data$clade_name,value=FALSE),]
  data_species <- data_species[grep('s__',data_species$clade_name,value=FALSE),]
  # data_species$clade_name <- sapply(strsplit(data_species$clade_name,'\\|'),'[',7)
  # data_species$clade_name <- sub("s__", "", data_species$clade_name)
  data_species <- as.data.frame(t(data_species))
  colnames(data_species) <- data_species['clade_name',]
  data_species <- data_species[-1,]
  
  if(length(unique(colnames(data_species)))<ncol(data_species)){
    message('dumplicated species!')
    species_unique <- unique(colnames(data_species))
    feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(data_species),ncol = length(species_unique)))
    colnames(feat_unique) <- species_unique
    rownames(feat_unique) <- rownames(data_species)
    for (species_t in species_unique) {
      match_t <- colnames(data_species)==species_t
      if(sum(match_t)>1){
        feat_t <- data_species[,colnames(data_species)==species_t]
        feat_unique[,species_t] <- rowSums(feat_t)
      }else{
        feat_unique[,species_t] <- data_species[,species_t]
      }
    }
    data_species <- feat_unique
  }
  # rowSums(feat_unique)
  # rowSums(feat)
  if(name_exclude==T){
    rownames(data_species) <- sapply(strsplit(rownames(data_species),'\\.'),'[',1)
  }
  data_species <- as.data.frame(t(data_species))
  # library(dplyr)
  data_species %>% mutate(across(where(is.character), as.numeric)) -> data_species
  return(data_species)
}

# SGB 
metaphlan4.species.name.SGB <- function(data,name_exclude=T){
  data_species <- data[grep('t__SGB',data$clade_name,value=FALSE),]
  data_species <- data_species[grep('k__Bacteria',data_species$clade_name,value=FALSE),]
  # data_species$clade_name <- sapply(strsplit(data_species$clade_name,'\\|'),'[',7)
  # data_species$clade_name <- sub("s__", "", data_species$clade_name)
  data_species <- as.data.frame(t(data_species))
  colnames(data_species) <- data_species['clade_name',]
  data_species <- data_species[-1,]
  
  if(length(unique(colnames(data_species)))<ncol(data_species)){
    message('dumplicated species!')
    species_unique <- unique(colnames(data_species))
    feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(data_species),ncol = length(species_unique)))
    colnames(feat_unique) <- species_unique
    rownames(feat_unique) <- rownames(data_species)
    for (species_t in species_unique) {
      match_t <- colnames(data_species)==species_t
      if(sum(match_t)>1){
        feat_t <- data_species[,colnames(data_species)==species_t]
        feat_unique[,species_t] <- rowSums(feat_t)
      }else{
        feat_unique[,species_t] <- data_species[,species_t]
      }
    }
    data_species <- feat_unique
  }
  # rowSums(feat_unique)
  # rowSums(feat)
  if(name_exclude==T){
    rownames(data_species) <- sapply(strsplit(rownames(data_species),'\\.'),'[',1)
  }
  data_species <- as.data.frame(t(data_species))
  # library(dplyr)
  data_species %>% mutate(across(where(is.character), as.numeric)) -> data_species
  return(data_species)
}
# beta.pcoa ----
sign.comp <- function(list) {
  pheno <- unique(list)
  result <- list()
  if(length(pheno)!=2){
    for (i in 1:(length(pheno) - 1)) {
      for (j in (i + 1):length(pheno)) {
        comp <- list(c(pheno[i], pheno[j]))
        result <- append(result, comp)
      }
    }
  }else{
    result <- append(result,list(c(pheno[1], pheno[2])))
  }
  
  return(result)
}
beta.pcoa <- function(dis_mat, group, groupID = "group",comp=NA, ellipse = T, label = F, PCo = 12, colors=NULL,p.all0=F ) {
  library(tidyverse)
  library(vegan)
  library(ape)
  library(phyloseq)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
  
  # group：一列sample，一列group
  dis_mat <- dis_mat[rownames(dis_mat) %in% group$sample,colnames(dis_mat) %in% group$sample]
  # PCoA
  pcoa <- cmdscale(dis_mat, k = 3, eig = T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  eig <- pcoa$eig
  points <- merge(points, group[, c("sample", groupID)], by.x = "row.names", by.y = "sample")
  rownames(points) <- points[, 1]
  points <- points[, -1]
  colnames(points) <- c("x", "y", "z", "group")
  if(is.na(comp) %>% unique){
    comp <- sign.comp(list = unique(group[,groupID]))
  }
  points$group <- factor(points$group, levels = group_all)
  # print(colors)
  # colors <- pal_lancet(palette = c("lanonc"), alpha = 0.8)(length(unique(group[,groupID])))
  if(is.null(colors)){
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    colors<-getPalette(length(unique(group[,groupID])))
  }
  
  
  
  # "#E41A1C" "#3A85A8" "#629363" "#C4625D" "#FFC81D" "#BF862B" "#EB7AA9" "#999999" "#D8D8D8"
  # colors<-c('#E41A1C','#4A72A6')
  # colors<-c('#E41A1C','#4A72A6')
  # colors<-c('#4A72A6','#3F7100','#E41A1C')
  
  # if (length(unique(points$group)) == 3) {
  #   colors <- c(colors[2], colors[3], colors[1])
  # }
  # if (length(unique(points$group)) == 2) {
  #   # colors <- c(colors[3], colors[2])
  #   colors <- c("#EEB043","#6B7BB3")
  #   
  # }
  
  
  
  # ADONIS
  group <- group[group$sample %in% rownames(dis_mat), ]
  anonis.result <- adonis2(dis_mat ~ group, data = group, permutations = 999)
  
  
  # 按1、2轴绘图
  if (PCo == 12) {
    p <- ggplot(points, aes(x = x, y = y, color = group)) +#, shape = group
      labs(x = "", y = "") +
      theme_bw() +
      guides(fill = "none") +
      geom_point(aes(color = group), alpha = .7, size = 2)#shape = group
    # labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
    #      y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按1、3轴绘图
  if (PCo == 13) {
    p <- ggplot(points, aes(x = x, y = z, color = group)) +
      labs(x = "", y = "") +
      theme_bw() +
      geom_point(aes(shape = group), alpha = .7, size = 2)
    # labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
    #      y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  # 按2、3轴绘图
  if (PCo == 23) {
    p <- ggplot(points, aes(x = y, y = z, color = group)) +
      labs(x = "", y = "") +
      theme_bw() +
      geom_point(aes(shape = group), alpha = .7, size = 2)
    # labs(x=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
    #      y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""), color=groupID)
  }
  
  p <- p +
    labs(title = "PCoA") +
    theme_bw() +
    geom_vline(xintercept = 0, lty = "dashed") +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) +
    scale_color_manual(name = "Group", values = colors) +
    scale_shape_discrete(name = "Group")
  text <- substitute(paste("PCoA ( Adonis ", R^2, " = ", m, ", p = ", n, ")"), list(m = round(anonis.result$R2[1], 2), n = anonis.result$`Pr(>F)`[1]))
  p <- p + labs(title = text)
  p
  # stat_ellipse(data=points,geom = "polygon",aes(fill=group),alpha=0.3)
  
  # geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")+
  # geom_vline(aes(xintercept=0),colour="#990000",linetype="dashed")
  
  p.x <- ggplot(data = points, aes(x = group, y = x, fill = group)) +
    geom_boxplot() +
    labs(y = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""), x = "") +
    scale_fill_manual(name = "Group", values = colors) +
    theme_bw() +
    coord_flip() +
    theme(text = element_text(size = 15), axis.ticks.y = element_blank()) +
    guides(fill = "none") +
    scale_x_discrete(labels = NULL) +
    # stat_compare_means(comparisons = lapply(comp,as.character), method = "wilcox.test")
    geom_signif(
      comparisons = lapply(comp,as.character),
      step_increase = 0.1,
      map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "."=0.1, "ns."=2),
      test = wilcox.test,
      test.args = c("two.sided"),
      tip_length = 0
    )
  
  p.y <- ggplot(data = points, aes(x = group, y = y, fill = group)) +
    geom_boxplot() +
    labs(y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = ""), x = "") +
    scale_fill_manual(name = "Group", values = colors) +
    theme_bw() +
    theme(text = element_text(size = 15)) +
    guides(fill = "none") +
    scale_x_discrete(labels = NULL) +
    geom_signif(
      comparisons = lapply(comp,as.character),
      step_increase = 0.1,
      map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "."=0.1, "ns."=2),
      test = wilcox.test, 
      test.args = c("two.sided"),
      tip_length = 0
    )
  p.z <- ggplot(data = points, aes(x = group, y = z, fill = group)) +
    geom_boxplot() +
    labs(y = paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits = 4), "%)", sep = ""), x = "") +
    scale_fill_manual(name = "Group", values = colors) +
    theme_bw() +
    theme(text = element_text(size = 15), axis.ticks.x = element_blank()) +
    guides(fill = "none") +
    scale_x_discrete(labels = NULL) +
    geom_signif(
      comparisons = lapply(comp,as.character),
      step_increase = 0.1,
      map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "."=0.1, "ns."=2),
      test = wilcox.test, 
      test.args = c("two.sided"),
      tip_length = 0
    )
  
  
  
  # 是否添加置信椭圆
  if (ellipse == T) {
    p <- p + stat_ellipse(level = 0.95)
  }
  # 是否显示样本标签
  if (label == T) {
    p <- p + geom_text_repel(label = paste(rownames(points)), size = 2.5)
  }
  
  if (PCo == 12) {
    p.y <- p.y + theme(axis.ticks.x = element_blank())
    p.all <- (p.y + p + plot_layout(nrow = 1, widths = c(1, 5))) - (plot_spacer() + p.x + plot_layout(nrow = 1, widths = c(1, 5))) + plot_layout(ncol = 1, heights = c(5, 1))
  }
  if (PCo == 13) {
    p.all <- (p.z + p + plot_layout(nrow = 1, widths = c(1, 5))) - (plot_spacer() + p.x + plot_layout(nrow = 1, widths = c(1, 5))) + plot_layout(ncol = 1, heights = c(5, 1))
  }
  if (PCo == 23) {
    p.y <- p.y + coord_flip() + theme(axis.ticks.y = element_blank())
    p.all <- (p.z + p + plot_layout(nrow = 1, widths = c(1, 5))) - (plot_spacer() + p.y + plot_layout(nrow = 1, widths = c(1, 5))) + plot_layout(ncol = 1, heights = c(5, 1))
  }
  
  p.all <- p.all + theme(text = element_text(family = "serif"))
  if(p.all0==T){
    return(p.all)
  }else{
    return(p)
  }
}

# marker ----
marker.identify.abundance <- function(feat,meta,level='genus',lda_cutoff=2){
  library(microbiomeMarker)
  library(phyloseq)
  library(ggplot2)
  library(RMySQL)
  
  otu_mat <- feat
  samples <- meta
  #筛选符合特征数大于2的样本
  #otu_mat <- otu_mat[rowSums(otu_mat>0)>2,]
  # 筛选特征
  otu_mat <- otu_mat[,colSums(otu_mat>0)>3]
  #同时改变meta
  samples <- samples[rownames(otu_mat),]
  #make phyloseq
  #使得otu_mat行名是1234... 列名是runid
  otu_mat <- t(otu_mat)
  tax_mat <- rownames(otu_mat)
  rownames(otu_mat) <- 1:length(tax_mat)
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  rownames(tax_mat) <- 1:nrow(tax_mat)
  if(level=='genus'){
    colnames(tax_mat) <- "Genus"
  }
  if(level=='species'){
    colnames(tax_mat) <- "Species"
  }
  
  ## combined to phyloseq format 
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(samples)
  phylodat <- phyloseq(OTU, TAX, samples)
  #lefse
  lefse.obj <- run_lefse(
    phylodat, 
    group = "Group", 
    multigrp_strat = TRUE,
    lda_cutoff = lda_cutoff #2021.8.12新增：不加此参数会默认为2
  )
  return(lefse.obj)
}
