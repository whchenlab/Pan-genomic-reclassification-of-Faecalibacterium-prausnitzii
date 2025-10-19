#### fun ----
source('./00_code/00_FUN/library_R_packages.R')
source('./00_code/00_FUN/siamcat_models.R')
source('./00_code/00_FUN/siamcat_models_adj.R')
source('./00_code/00_FUN/my_lodo.R')
source('./00_code/00_FUN/my_lodo.adj.R')
source('./00_code/00_FUN/00_Fun.R')
source('./00_code/00_FUN/CCM & SCM Fun.R')
setwd('E:/OneDrive - hust.edu.cn/项目/chenlab/2025-FP-FD/')



# 柱状图 F genus -----
type='gut'
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggplotify)
library(grid)

if(type!='gut'){
  sum <- read.csv("01_data/cluster_number.csv")
  
  # sum_df <- data.frame(sum_62=c(95,95,85),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"))))
  sum_df <- data.frame(sum_62=sum$Number,Group = sum$Type.Strain)
  sum_df$sum_62 <- as.numeric(sum_df$sum_62)
  
  sum_df$Group <- sapply(strsplit(sum_df$Group, "_"), function(x) paste0(x[1]))
  sum_df$Group[sum_df$Group=="<ce><de>"] <- "Others"
  # sum(sum_df[sum_df$Group=="Others",]$sum_62)
  sum_df <- sum_df[sum_df$Group!="Others",]
  sum_df$Group[sum_df$Group=="F.taiwanense"] <- "F.taiwanense; F.faecis"
  
  
  sum_df$Group <- factor(sum_df$Group, levels = sum_df$Group)
}else{
  sum <- read.table("01_data/SpeciesNumber_Genus.txt")
  colnames(sum) <- c('Type.strain','Number')
  sum <- rbind(sum,c('Others',sum(sum$Number[sum$Type.strain=='Unknown'])))
  sum <- sum[-which(sum$Type.strain=='Unknown'),]
  sum_df <- data.frame(sum_62=sum$Number,Group = sum$Type.strain)
  sum_df$sum_62 <- as.numeric(sum_df$sum_62)
  
  sum_df$Group <- factor(sum_df$Group, levels = sum_df$Group)
}

# 假设sum_df有Group和sum_62两列数据
ggplot(sum_df, aes(x = Group, y = sum_62, fill = sum_62)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic2() +
  labs(title = "", x = "", y = "Genome number", fill = " ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 添加渐变色填充
  # scale_fill_gradient(
  #   low = "#66C2A4",    # 渐变色起始颜色（低数值）
  #   high = "#00441B",   # 渐变色结束颜色（高数值）
  #   guide = guide_colorbar(barwidth = 1.5, barheight = 10)
  # ) +
  
  scale_fill_gradient(
    low = "#CCEBC5",    # 渐变色起始颜色（低数值）
    high = "#084081",   # 渐变色结束颜色（高数值）
    guide = guide_colorbar(barwidth = 1.5, barheight = 10)
  ) +
  # 添加柱状图上方的数字标注
  geom_text(
    aes(label = sum_62),
    vjust = -0.5,       # 调整文本位置在柱子上方
    size = 4,           # 文本大小
    color = "black"     # 文本颜色
  )    
#标注最高
# geom_hline(yintercept = max(sum_62), color = "red") +
# annotate("text", x = 1, y = max(sum_62) + 0.5, label = paste("   ", max(sum_62)), color = "red")
ggsave(filename = "ani_genome_number_genus_gut.pdf", width = 5, height = 5)



# 柱状图 Fp species  -----
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggplotify)
library(grid)
# sum <- read.csv("01_data/FP.csv")
sum <- read.table("01_data/SpeciesNumber_FP.txt")
colnames(sum) <- c('Type.strain','Number')
# sum$Type <- NULL
# sum <- rbind(sum[grep('F',sum$Type.strain,invert = F),], data.frame( Number=sum(sum$Number[grep('F',sum$Type.strain,invert = T)]),Type.strain="Others"))
sum <- rbind(sum,c('Others',sum(sum$Number[sum$Type.strain=='Unknown'])))
sum <- sum[-which(sum$Type.strain=='Unknown'),]

sum_df <- data.frame(sum_62=sum$Number,Group = sum$Type.strain)
sum_df$sum_62 <- as.numeric(sum_df$sum_62)

sum_df$Group <- factor(sum_df$Group, levels = sum_df$Group)
# 假设sum_df有Group和sum_62两列数据
ggplot(sum_df, aes(x = Group, y = sum_62, fill = sum_62)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic2() +
  labs(title = "", x = "", y = "Genome number", fill = " ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face = 'italic')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 添加渐变色填充
  # scale_fill_gradient(
  #   low = "#66C2A4",    # 渐变色起始颜色（低数值）
  #   high = "#00441B",   # 渐变色结束颜色（高数值）
  #   guide = guide_colorbar(barwidth = 1.5, barheight = 10)
  # ) +
  
  scale_fill_gradient(
    low = "#CCEBC5",    # 渐变色起始颜色（低数值）
    high = "#084081",   # 渐变色结束颜色（高数值）
    guide = guide_colorbar(barwidth = 1.5, barheight = 10)
  ) +
  # 添加柱状图上方的数字标注
  geom_text(
    aes(label = sum_62),
    vjust = -0.5,       # 调整文本位置在柱子上方
    size = 4,           # 文本大小
    color = "black"     # 文本颜色
  )    
#标注最高
# geom_hline(yintercept = max(sum_62), color = "red") +
# annotate("text", x = 1, y = max(sum_62) + 0.5, label = paste("   ", max(sum_62)), color = "red")
ggsave(filename = "FP_ani_genome_number_gut.pdf", width = 5, height = 5)


# 堆积柱状图 FAAH.bar ----
sum_df <- data.frame(sum_62=c(94,10,82,0,82,0),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"), 2)),
                     Category=paste0(rep(c("No","Yes"), each=3)))
sum_df$sum_62 <- as.numeric(sum_df$sum_62)
sum_df$Group <- factor(sum_df$Group, levels = c("F.prausnitzii","F.longum","F.duncaniae"))
ggplot(sum_df, aes(x = Group, y = sum_62, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  labs(title = "", x = "", y = "Genome number", fill = "FAAH") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 标注每组的总和
  geom_text(
    aes(label = sum_62, group = Category),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "white"
  ) +
  # 标注所有分组中的最大值
  # geom_hline(yintercept = max(tapply(sum_df$sum_62, sum_df$Group, sum)), color = "red") +
  # annotate(
  #   "text",
  #   x = 1,
  #   y = max(tapply(sum_df$sum_62, sum_df$Group, sum)) + 0.5,
  #   label = paste("Max: ", max(tapply(sum_df$sum_62, sum_df$Group, sum))),
  #   color = "red"
  # )+
  # 设置自定义颜色方案
  scale_fill_manual(
    values = c(
      "#FDBF6F", "#FF7F00"
    )
  )+
  # scale_fill_manual(
  #   values = c(
  #     "#56B4E9", "#009E73", "#F0E442",
  #     "#D55E00", "#CC79A7", "#E69F00",
  #     "#0072B2", "#999999"
  #   )
  # )+
  # scale_y_continuous(breaks = seq(4.5,6,0.5))+
  theme_bw()+
  theme(legend.position = 'right')+
  theme(
    text = element_text(family = "",size = 12,face = 'plain'),
    # axis.line=element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text( hjust=0.5, vjust=0,family = '',size = 10,face = 'italic'),#angle=90,
    # axis.text.x =element_blank(),
    # axis.ticks.x = element_blank(),
    # panel.grid=element_blank(),
    # panel.border=element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 12),
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),     legend.position = 'right'
  )
ggsave(filename = "FAAH.bar.pdf", width = 3.8, height = 3.5)



# data_df <- data.frame(FP=c(95,94,95,95,94,0,1,0,0,1),
#            FL=c(90,95,95,95,95,5,0,0,0,0),
#            FD=c(80,85,83,81,83,5,0,2,4,2),
#            gene=rep(c("bcd","but","cro","hbd","thl"),2),
#            Category=rep(c("Yes","No"),each=5))

# 堆积柱状图 FAAH.bar others ----
#### 读取单个 sheet
library(readxl)
sum_df_new <- read_excel("./01_data/FAAH.gut.all.xlsx",sheet = 1)
sum_df_new <-  sum_df_new[,-4]
colnames(sum_df_new) <- c('Group','Yes','No')
sum_df <- sum_df_new %>% gather(Category,sum_62,-'Group')
sum_df_new$Group[sum_df_new$Group=='proteins_fb'] <- 'F.butyricigenerans'
sum_df_new$Group[sum_df_new$Group=='proteins_ff'] <- 'F.faecis&F.taiwanense'
sum_df_new$Group[sum_df_new$Group=='proteins_fh'] <- 'F.hattorii'
sum_df_new$Group[sum_df_new$Group=='proteins_fi'] <- 'F.intestinale'
sum_df_new$Group[sum_df_new$Group=='proteins_fw'] <- 'F.wellingii'
sum_df_new <- rbind(sum_df_new,c('Unknown',sum(sum_df_new$Yes[grep('proteins_un',sum_df_new$Group)]),
                                 sum(sum_df_new$No[grep('proteins_un',sum_df_new$Group)])
))
sum_df_new$Yes <- as.numeric(sum_df_new$Yes)
sum_df_new$No <- as.numeric(sum_df_new$No)
sum_df_new <- sum_df_new[-grep('proteins_un',sum_df_new$Group),]
sum_df_new$sum <- sum_df_new$No+sum_df_new$Yes
sum_df_new <- arrange(sum_df_new,desc(sum))
sum_df_new$Group <- factor(sum_df_new$Group,levels = sum_df_new$Group)
sum_df_new$sum <- NULL
# sum_df <- data.frame(sum_62=c(95,10,85,0,85,0),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"), 2)),
#                      Category=paste0(rep(c("No","Yes"), each=3)))
# sum_df$sum_62 <- as.numeric(sum_df$sum_62)
# sum_df$Group <- factor(sum_df$Group, levels = c("F.prausnitzii","F.longum","F.duncaniae"))
sum_df <- gather(sum_df_new,'Category','sum_62',-"Group")
ggplot(sum_df, aes(x = Group, y = sum_62, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  labs(title = "", x = "", y = "Genome number", fill = "FAAH") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 标注每组的总和
  geom_text(
    aes(label = sum_62, group = Category),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "white"
  ) +
  # 标注所有分组中的最大值
  # geom_hline(yintercept = max(tapply(sum_df$sum_62, sum_df$Group, sum)), color = "red") +
  # annotate(
  #   "text",
  #   x = 1,
  #   y = max(tapply(sum_df$sum_62, sum_df$Group, sum)) + 0.5,
  #   label = paste("Max: ", max(tapply(sum_df$sum_62, sum_df$Group, sum))),
  #   color = "red"
  # )+
  # 设置自定义颜色方案
  scale_fill_manual(
    values = c(
      "#FDBF6F", "#FF7F00"
    )
  )+
  # scale_fill_manual(
  #   values = c(
  #     "#56B4E9", "#009E73", "#F0E442",
  #     "#D55E00", "#CC79A7", "#E69F00",
  #     "#0072B2", "#999999"
  #   )
  # )+
  # scale_y_continuous(breaks = seq(4.5,6,0.5))+
  theme_bw()+
  theme(legend.position = 'right')+
  theme(
    text = element_text(family = "",size = 12,face = 'plain'),
    # axis.line=element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle=90, hjust=0.5, vjust=0,family = '',size = 10,face = 'italic'),#angle=90,
    # axis.text.x =element_blank(),
    # axis.ticks.x = element_blank(),
    # panel.grid=element_blank(),
    # panel.border=element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 12),
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),     legend.position = 'right'
  )

ggsave(filename = "FAAH.bar.all.pdf", width = 5, height = 5)



# data_df <- data.frame(FP=c(95,94,95,95,94,0,1,0,0,1),
#            FL=c(90,95,95,95,95,5,0,0,0,0),
#            FD=c(80,85,83,81,83,5,0,2,4,2),
#            gene=rep(c("bcd","but","cro","hbd","thl"),2),
#            Category=rep(c("Yes","No"),each=5))



# 堆积柱状图 but ----
library(ggplot2)
library(dplyr)
library(tidyr)

library(readxl)
sum_df_new <- read_excel("./01_data/Butyrate.xlsx",sheet = 1)
sum_df_new <-  sum_df_new[,-5]
colnames(sum_df_new) <- c('Gene','Group','Yes','No')
sum_df_new$Gene <- gsub('.faa','',sum_df_new$Gene)
sum_df_new$Gene <- factor(sum_df_new$Gene, levels = c("thl", "hbd", "cro", "bcd", "but"))


# sum_df_new$Group[sum_df_new$Group=='proteins_fb'] <- 'F.butyricigenerans'
# sum_df_new$Group[sum_df_new$Group=='proteins_ff'] <- 'F.faecis&F.taiwanense'
# sum_df_new$Group[sum_df_new$Group=='proteins_fh'] <- 'F.hattorii'
# sum_df_new$Group[sum_df_new$Group=='proteins_fi'] <- 'F.intestinale'
# sum_df_new$Group[sum_df_new$Group=='proteins_fw'] <- 'F.wellingii'
sum_df_new$Group[sum_df_new$Group=='proteins_fp'] <- 'F.prausnitzii'
sum_df_new$Group[sum_df_new$Group=='proteins_fl'] <- 'F.longum'
sum_df_new$Group[sum_df_new$Group=='proteins_fd'] <- 'F.duncaniae'


# sum_df_new <- rbind(sum_df_new,sum_df_new %>% 
#                       filter(grepl('proteins_un',Group))  %>% 
#                       group_by(Gene) %>% 
#                       summarise(Group='Unknown',Yes= sum(Yes), No = sum(No)))
# sum_df_new <- rbind(sum_df_new,c('Unknown',sum(sum_df_new$Yes[grep('proteins_un',sum_df_new$Group)]),
#                                  sum(sum_df_new$No[grep('proteins_un',sum_df_new$Group)])
# ))
sum_df_new$Yes <- as.numeric(sum_df_new$Yes)
sum_df_new$No <- as.numeric(sum_df_new$No)
# sum_df_new <- sum_df_new[-grep('proteins_un',sum_df_new$Group),]
# sum_df_new$sum <- sum_df_new$No+sum_df_new$Yes
# sum_df_new <- arrange(sum_df_new,sum)
sum_df_new$Group <- factor(sum_df_new$Group,levels = c('F.duncaniae','F.longum','F.prausnitzii'))
sum_df_new$sum <- NULL
# sum_df <- data.frame(sum_62=c(95,10,85,0,85,0),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"), 2)),
#                      Category=paste0(rep(c("No","Yes"), each=3)))
# sum_df$sum_62 <- as.numeric(sum_df$sum_62)
# sum_df$Group <- factor(sum_df$Group, levels = c("F.prausnitzii","F.longum","F.duncaniae"))
sum_df <- gather(sum_df_new,'Category','sum_62',-c("Gene","Group"))
# 创建堆积柱状图
sum_df$label <- ifelse(sum_df$sum_62!=0,sum_df$sum_62,'')

# 设置自定义颜色
custom_colors <- c("#CAB2D6", "#6A3D9A")
# custom_colors <- c("#6A3D9A", "#CAB2D6")

stacked_barplot <- ggplot(sum_df, aes(x = sum_62, y =Group , fill = Category)) +
  geom_bar(stat = "identity", width = 0.25) +
  # facet_wrap(~ Bacteria_Full, nrow = 1) +
  facet_grid( ~ Gene) +
  theme_minimal() +
  labs(x = "", y = "Genome number", fill = "Category") +
  theme(
    text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10,face = "italic"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 标注每组的数值
  geom_text(
    aes(label = label, group = Category),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white"
  ) +
  # 设置自定义颜色方案
  scale_fill_manual(values = custom_colors)

# 显示图形
print(stacked_barplot)

# ggsave(filename = "Butyrate.bar.pdf", width = 6, height = 2.7)
ggsave(filename = "Butyrate.bar_deep_gut.pdf", width = 9, height = 2.3)




# # 原始数据
# data_df <- data.frame(
#   FP = c(95, 94, 95, 95, 94, 0, 1, 0, 0, 1),
#   FL = c(90, 95, 95, 95, 95, 5, 0, 0, 0, 0),
#   FD = c(80, 85, 83, 81, 83, 5, 0, 2, 4, 2),
#   gene = rep(c("bcd", "but", "cro", "hbd", "thl"), 2),
#   Category = rep(c("Yes", "No"), each = 5)
# )
# 
# data_df$gene <- factor(data_df$gene, levels = c("thl", "hbd", "cro", "bcd", "but"))
# 
# # 转换数据格式为长格式
# data_long <- data_df %>%
#   pivot_longer(cols = c(FP, FL, FD), names_to = "Bacteria", values_to = "Count")
# 
# # 计算每个基因-细菌组合的总和
# total_counts <- data_long %>%
#   group_by(gene, Bacteria) %>%
#   summarise(Total = sum(Count))
# 
# # 合并原始数据和总计数
# plot_data <- data_long %>%
#   left_join(total_counts, by = c("gene", "Bacteria")) %>%
#   mutate(Percentage = Count / Total * 100)
# 
# # 映射Bacteria缩写到全称
# bacteria_full_names <- c(
#   "FP" = "F.prausnitzii",
#   "FL" = "F.longum",
#   "FD" = "F.duncaniae"
# )
# plot_data$Bacteria_Full <- bacteria_full_names[plot_data$Bacteria]
# plot_data$Bacteria_Full <- factor(plot_data$Bacteria_Full, 
#                                   levels = c("F.duncaniae", "F.longum", "F.prausnitzii"))
# 
# # 设置自定义颜色
# custom_colors <- c("#CAB2D6", "#6A3D9A")
# custom_colors <- c("#6A3D9A", "#CAB2D6")
# 
# # 创建堆积柱状图
# stacked_barplot <- ggplot(plot_data, aes(x = Count, y =Bacteria_Full , fill = Category)) +
#   geom_bar(stat = "identity", width = 0.2) +
#   # facet_wrap(~ Bacteria_Full, nrow = 1) +
#   facet_grid( ~ gene) +
#   theme_minimal() +
#   labs(x = "", y = "Genome number", fill = "Category") +
#   theme(
#     text = element_text(size = 12),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 10,face = "italic"),
#     strip.text = element_text(size = 12, face = "bold"),
#     legend.position = "right",
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) +
#   # 标注每组的数值
#   geom_text(
#     aes(label = Count, group = Category),
#     position = position_stack(vjust = 0.5),
#     size = 4,
#     color = "white"
#   ) +
#   # 设置自定义颜色方案
#   scale_fill_manual(values = custom_colors)
# 
# # 显示图形
# print(stacked_barplot)

ggsave(filename = "Butyrate.bar.pdf", width = 6, height = 2.7)
ggsave(filename = "Butyrate.bar_deep.pdf", width = 9, height = 2.7)



# 堆积柱状图 but others ----
#### 读取单个 sheet
library(readxl)
sum_df_new <- read_excel("./01_data/Butyrate.gut.all.xlsx",sheet = 1)
sum_df_new <-  sum_df_new[,-5]
colnames(sum_df_new) <- c('Gene','Group','Yes','No')
sum_df_new$Gene <- gsub('.faa','',sum_df_new$Gene)
sum_df_new$Gene <- factor(sum_df_new$Gene, levels = c("thl", "hbd", "cro", "bcd", "but"))


sum_df_new$Group[sum_df_new$Group=='proteins_fb'] <- 'F.butyricigenerans'
sum_df_new$Group[sum_df_new$Group=='proteins_ff'] <- 'F.faecis&F.taiwanense'
sum_df_new$Group[sum_df_new$Group=='proteins_fh'] <- 'F.hattorii'
sum_df_new$Group[sum_df_new$Group=='proteins_fi'] <- 'F.intestinale'
sum_df_new$Group[sum_df_new$Group=='proteins_fw'] <- 'F.wellingii'
sum_df_new <- rbind(sum_df_new,sum_df_new %>% 
                      filter(grepl('proteins_un',Group))  %>% 
                      group_by(Gene) %>% 
                      summarise(Group='Unknown',Yes= sum(Yes), No = sum(No)))
# sum_df_new <- rbind(sum_df_new,c('Unknown',sum(sum_df_new$Yes[grep('proteins_un',sum_df_new$Group)]),
#                                  sum(sum_df_new$No[grep('proteins_un',sum_df_new$Group)])
# ))
sum_df_new$Yes <- as.numeric(sum_df_new$Yes)
sum_df_new$No <- as.numeric(sum_df_new$No)
sum_df_new <- sum_df_new[-grep('proteins_un',sum_df_new$Group),]
sum_df_new$sum <- sum_df_new$No+sum_df_new$Yes
sum_df_new <- arrange(sum_df_new,sum)
sum_df_new$Group <- factor(sum_df_new$Group,levels = unique(sum_df_new$Group))
sum_df_new$sum <- NULL
# sum_df <- data.frame(sum_62=c(95,10,85,0,85,0),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"), 2)),
#                      Category=paste0(rep(c("No","Yes"), each=3)))
# sum_df$sum_62 <- as.numeric(sum_df$sum_62)
# sum_df$Group <- factor(sum_df$Group, levels = c("F.prausnitzii","F.longum","F.duncaniae"))
sum_df <- gather(sum_df_new,'Category','sum_62',-c("Gene","Group"))

# 设置自定义颜色
custom_colors <- c("#CAB2D6", "#6A3D9A")
# custom_colors <- c("#6A3D9A", "#CAB2D6")

# 创建堆积柱状图
sum_df$label <- ifelse(sum_df$sum_62!=0,sum_df$sum_62,'')
stacked_barplot <- ggplot(sum_df, aes(x = sum_62, y =Group , fill = Category)) +
  geom_bar(stat = "identity", width = 0.3) +
  # facet_wrap(~ Bacteria_Full, nrow = 1) +
  facet_grid( ~ Gene) +
  theme_classic() +
  labs(x = "", y = "Genome number", fill = "Category") +
  theme(
    text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10,face = "italic"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 标注每组的数值
  geom_text(
    aes(label = label, group = Category),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white"
  ) +
  # 设置自定义颜色方案
  scale_fill_manual(values = custom_colors)

# 显示图形
print(stacked_barplot)

# ggsave(filename = "Butyrate.bar.pdf", width = 6, height = 2.7)
ggsave(filename = "Butyrate.bar_deep_all.pdf", width = 12, height = 3.5)

# GalNAc N-乙酰半乳糖胺 bar plot ----

sum_df <- data.frame(sum_62=c(86,10,75,8,82,7),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"), 2)),
                     Category=paste0(rep(c("No","Yes"), each=3)))
sum_df$sum_62 <- as.numeric(sum_df$sum_62)
sum_df$Group <- factor(sum_df$Group, levels = c("F.prausnitzii","F.longum","F.duncaniae"))
ggplot(sum_df, aes(x = Group, y = sum_62, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  labs(title = "", x = "", y = "Genome number", fill = "GalNAc") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 标注每组的总和
  geom_text(
    aes(label = sum_62, group = Category),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "white"
  ) +
  # 标注所有分组中的最大值
  # geom_hline(yintercept = max(tapply(sum_df$sum_62, sum_df$Group, sum)), color = "red") +
  # annotate(
  #   "text",
  #   x = 1,
  #   y = max(tapply(sum_df$sum_62, sum_df$Group, sum)) + 0.5,
  #   label = paste("Max: ", max(tapply(sum_df$sum_62, sum_df$Group, sum))),
  #   color = "red"
  # )+
  # 设置自定义颜色方案
  scale_fill_manual(
    values = c(
      # "#FDBF6F", "#FF7F00"
      # "#A6CEE3", "#1F78B4"
      # "#B2DF8A" ,"#33A02C"
      "#FB9A99" ,"#E31A1C"
    )
  )+
  # scale_fill_manual(
  #   values = c(
  #     "#56B4E9", "#009E73", "#F0E442",
  #     "#D55E00", "#CC79A7", "#E69F00",
  #     "#0072B2", "#999999"
  #   )
  # )+
  # scale_y_continuous(breaks = seq(4.5,6,0.5))+
  theme_bw()+
  theme(legend.position = 'right')+
  theme(
    text = element_text(family = "",size = 12,face = 'plain'),
    # axis.line=element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text( hjust=0.5, vjust=0,family = '',size = 10,face = 'italic'),#angle=90,
    # axis.text.x =element_blank(),
    # axis.ticks.x = element_blank(),
    # panel.grid=element_blank(),
    # panel.border=element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 12),
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),     legend.position = 'right'
  )
ggsave(filename = "GalNAc.bar.pdf", width = 3.8, height = 3.5)


# GalNAc N-乙酰半乳糖胺 bar others ----
#### 读取单个sheet
library(readxl)
sum_df_new <- read_excel("./01_data/dsv.gut.all.xlsx",sheet = 1)
sum_df_new <-  sum_df_new[,-4]
colnames(sum_df_new) <- c('Group','Yes','No')
sum_df <- sum_df_new %>% gather(Category,sum_62,-'Group')
sum_df_new$Group[sum_df_new$Group=='proteins_fb'] <- 'F.butyricigenerans'
sum_df_new$Group[sum_df_new$Group=='proteins_ff'] <- 'F.faecis&F.taiwanense'
sum_df_new$Group[sum_df_new$Group=='proteins_fh'] <- 'F.hattorii'
sum_df_new$Group[sum_df_new$Group=='proteins_fi'] <- 'F.intestinale'
sum_df_new$Group[sum_df_new$Group=='proteins_fw'] <- 'F.wellingii'
sum_df_new <- rbind(sum_df_new,c('Unknown',sum(sum_df_new$Yes[grep('proteins_un',sum_df_new$Group)]),
                                 sum(sum_df_new$No[grep('proteins_un',sum_df_new$Group)])
))
sum_df_new$Yes <- as.numeric(sum_df_new$Yes)
sum_df_new$No <- as.numeric(sum_df_new$No)
sum_df_new <- sum_df_new[-grep('proteins_un',sum_df_new$Group),]
sum_df_new <- sum_df_new[-grep('proteins_fp',sum_df_new$Group),]
sum_df_new <- sum_df_new[-grep('proteins_fl',sum_df_new$Group),]
sum_df_new <- sum_df_new[-grep('proteins_fd',sum_df_new$Group),]
sum_df_new$sum <- sum_df_new$No+sum_df_new$Yes
sum_df_new <- arrange(sum_df_new,desc(sum))
sum_df_new$Group <- factor(sum_df_new$Group,levels = sum_df_new$Group)
sum_df_new$sum <- NULL
# sum_df <- data.frame(sum_62=c(95,10,85,0,85,0),Group = paste0(rep(c("F.prausnitzii","F.longum","F.duncaniae"), 2)),
#                      Category=paste0(rep(c("No","Yes"), each=3)))
# sum_df$sum_62 <- as.numeric(sum_df$sum_62)
# sum_df$Group <- factor(sum_df$Group, levels = c("F.prausnitzii","F.longum","F.duncaniae"))
sum_df <- gather(sum_df_new,'Category','sum_62',-"Group")
ggplot(sum_df, aes(x = Group, y = sum_62, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  labs(title = "", x = "", y = "Genome number", fill = "GalNAc") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 标注每组的总和
  geom_text(
    aes(label = sum_62, group = Category),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "white"
  ) +
  # 标注所有分组中的最大值
  # geom_hline(yintercept = max(tapply(sum_df$sum_62, sum_df$Group, sum)), color = "red") +
  # annotate(
  #   "text",
  #   x = 1,
  #   y = max(tapply(sum_df$sum_62, sum_df$Group, sum)) + 0.5,
  #   label = paste("Max: ", max(tapply(sum_df$sum_62, sum_df$Group, sum))),
  #   color = "red"
  # )+
  # 设置自定义颜色方案
  scale_fill_manual(
    values = c(
      "#FB9A99" ,"#E31A1C"
    )
  )+
  # scale_fill_manual(
  #   values = c(
  #     "#56B4E9", "#009E73", "#F0E442",
  #     "#D55E00", "#CC79A7", "#E69F00",
  #     "#0072B2", "#999999"
  #   )
  # )+
  # scale_y_continuous(breaks = seq(4.5,6,0.5))+
  theme_bw()+
  theme(legend.position = 'right')+
  theme(
    text = element_text(family = "",size = 12,face = 'plain'),
    # axis.line=element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle=90, hjust=0.5, vjust=0,family = '',size = 10,face = 'italic'),#angle=90,
    # axis.text.x =element_blank(),
    # axis.ticks.x = element_blank(),
    # panel.grid=element_blank(),
    # panel.border=element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 12),
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),     legend.position = 'right'
  )

ggsave(filename = "GalNAc.bar.all.pdf", width = 5, height = 5)



# data_df <- data.frame(FP=c(95,94,95,95,94,0,1,0,0,1),
#            FL=c(90,95,95,95,95,5,0,0,0,0),
#            FD=c(80,85,83,81,83,5,0,2,4,2),
#            gene=rep(c("bcd","but","cro","hbd","thl"),2),
#            Category=rep(c("Yes","No"),each=5))

# 分类 热图展示  plot -----

# 1) NCBI ----
AUC_mat <- data.frame(Fp=c(111,0,0,0),FL=c(0,88,0,5),FD=c(0,0,81,0),Others=c(1,3,0,181))
rownames(AUC_mat) <- c("F.prausnitzii","F.longum","F.duncaniae","Others")
colnames(AUC_mat) <- c("F.prausnitzii","F.longum","F.duncaniae","Others")
datas <- AUC_mat
datas$TR <- rownames(datas)
datas$method <- 'HMM'
#col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
# col.scheme.heatmap= c('gray95','steelblue1','midnightblue')
col.scheme.heatmap= c('#F7FBFF','steelblue1','#08306B')
# col.scheme.heatmap= c('white','pink','red')

proj=unique(as.character(datas$TR))

data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
data_df=tibble(data_df)

g <- data_df %>% 
  mutate(study.train=TR) %>%
  filter(method == 'HMM') %>%
  mutate(study.test=factor(study.test, levels=proj)) %>% 
  mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0, 181)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        axis.text.x = element_text(size=10, face='italic'),
        axis.text.y = element_text(size=10, face='italic'),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  labs(fill='Number') +
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
g
ggsave(filename = "species_heatmap_ncbi.pdf", width = 5, height = 3.8)


# 2) BGI allfna2 ----

AUC_mat <- data.frame(Fp=c(47,0,0,0),FL=c(0,17,2,0),FD=c(0,0,37,0),Others=c(0,1,0,133))
rownames(AUC_mat) <- c("F.prausnitzii","F.longum","F.duncaniae","Others")
colnames(AUC_mat) <- c("F.prausnitzii","F.longum","F.duncaniae","Others")
datas <- AUC_mat
datas$TR <- rownames(datas)
datas$method <- 'HMM'
#col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
# col.scheme.heatmap= c('gray95','steelblue1','midnightblue')
col.scheme.heatmap= c('#F7FBFF','steelblue1','#08306B')
# col.scheme.heatmap= c('white','pink','red')

proj=unique(as.character(datas$TR))

data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
data_df=tibble(data_df)

g <- data_df %>% 
  mutate(study.train=TR) %>%
  filter(method == 'HMM') %>%
  mutate(study.test=factor(study.test, levels=proj)) %>% 
  mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0, 133)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        axis.text.x = element_text(size=10, face='italic'),
        axis.text.y = element_text(size=10, face='italic'),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  labs(fill='Number') +
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
g
ggsave(filename = "species_heatmap_BGI.pdf", width = 5, height = 3.8)

# 3) genus ----

AUC_mat <- data.frame(Faecalibacterium=c(11,0),non_Faecalibacterium=c(0,1740))
rownames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
colnames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
datas <- AUC_mat
datas$TR <- rownames(datas)
datas$method <- 'HMM'
#col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
# col.scheme.heatmap= c('gray95','steelblue1','midnightblue')
col.scheme.heatmap= c('#F7FBFF','steelblue1','#08306B')
# col.scheme.heatmap= c('white','pink','red')

proj=unique(as.character(datas$TR))

data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
data_df=tibble(data_df)

g1 <- data_df %>% 
  mutate(study.train=TR) %>%
  filter(method == 'HMM') %>%
  mutate(study.test=factor(study.test, levels=proj)) %>% 
  mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0, 1740)) + #max(data_df$AUC)
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        axis.text.x = element_text(size=10, face='italic'),
        axis.text.y = element_text(size=10, face='italic'),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  labs(fill='Number') +
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
# g
# ggsave(filename = "species_heatmap_ncbi_genus.pdf", width = 5, height = 3.8)


AUC_mat <- data.frame(Faecalibacterium=c(237,0),non_Faecalibacterium=c(0,0))
rownames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
colnames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
datas <- AUC_mat
datas$TR <- rownames(datas)
datas$method <- 'HMM'
#col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
# col.scheme.heatmap= c('gray95','steelblue1','midnightblue')
col.scheme.heatmap= c('#F7FBFF','steelblue1','#08306B')
# col.scheme.heatmap= c('white','pink','red')

proj=unique(as.character(datas$TR))

data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
data_df=tibble(data_df)

g2 <- data_df %>% 
  mutate(study.train=TR) %>%
  filter(method == 'HMM') %>%
  mutate(study.test=factor(study.test, levels=proj)) %>% 
  mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0, 1740)) + #max(data_df$AUC)
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        axis.text.x = element_text(size=10, face='italic'),
        axis.text.y = element_text(size=10, face='italic'),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  labs(fill='Number') +
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
# g
# ggsave(filename = "species_heatmap_ncbi_genus.pdf", width = 5, height = 3.8)

AUC_mat <- data.frame(Faecalibacterium=c(470,0),non_Faecalibacterium=c(0,0))
rownames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
colnames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
datas <- AUC_mat
datas$TR <- rownames(datas)
datas$method <- 'HMM'
#col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
# col.scheme.heatmap= c('gray95','steelblue1','midnightblue')
col.scheme.heatmap= c('#F7FBFF','steelblue1','#08306B')
# col.scheme.heatmap= c('white','pink','red')

proj=unique(as.character(datas$TR))

data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
data_df=tibble(data_df)

g3 <- data_df %>% 
  mutate(study.train=TR) %>%
  filter(method == 'HMM') %>%
  mutate(study.test=factor(study.test, levels=proj)) %>% 
  mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0, 1740)) + #max(data_df$AUC)
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        axis.text.x = element_text(size=10, face='italic'),
        axis.text.y = element_text(size=10, face='italic'),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  labs(fill='Number') +
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
# g
# ggsave(filename = "species_heatmap_ncbi_genus.pdf", width = 5, height = 3.8)

AUC_mat <- data.frame(Faecalibacterium=c(11,0),non_Faecalibacterium=c(0,329))
rownames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
colnames(AUC_mat) <- c("Faecalibacterium","non-Faecalibacterium")
datas <- AUC_mat
datas$TR <- rownames(datas)
datas$method <- 'HMM'
#col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
# col.scheme.heatmap= c('gray95','steelblue1','midnightblue')
col.scheme.heatmap= c('#F7FBFF','steelblue1','#08306B')
# col.scheme.heatmap= c('white','pink','red')

proj=unique(as.character(datas$TR))

data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
data_df=tibble(data_df)

g4 <- data_df %>% 
  mutate(study.train=TR) %>%
  filter(method == 'HMM') %>%
  mutate(study.test=factor(study.test, levels=proj)) %>% 
  mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0, 1740)) + #max(data_df$AUC)
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        axis.text.x = element_text(size=10, face='italic'),
        axis.text.y = element_text(size=10, face='italic'),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  labs(fill='Number') +
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
# g
# ggsave(filename = "species_heatmap_ncbi_genus.pdf", width = 5, height = 3.8)



library(ggarrange)
G <- ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2,
               labels = c("CGR2", "BGI", "NCBI", "HiBC"),
               common.legend = TRUE, legend = "right")
print(G)
ggsave(filename = "species_heatmap_genus.pdf", plot = G, width = 8, height = 6)





# ani matrix 聚类 heatmap----
ani <- read.table("01_data/ani_matrix_sgb.txt",  sep = "\t")
#长-》宽
library(tidyr)
ani_t <- ani[,-c(4,5)]
ani_long <- ani_t %>%
  spread(key = V1, value = V3)  
row.names(ani_long) <- ani_long$V2
ani_long$V2 <- NULL
rownames(ani_long) <- gsub(".fna", "", rownames(ani_long))
colnames(ani_long) <- gsub(".fna", "", colnames(ani_long))

# ani_matrix_typestrain_SGB.pdf
ani_long_t <- ani_long[grep("F.", row.names(ani_long)), grep("SGB", colnames(ani_long))] 
rownames(ani_long_t) <- sapply(strsplit(rownames(ani_long_t), "_"), function(x) paste(x[1], collapse = "_"))
library(pheatmap)
library(RColorBrewer)
heatmap_colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
#只显示大于95的值
ani_long_t_t <- ani_long_t
ani_long_t_t <- apply(ani_long_t_t, 2, function(x) round(x, 1))
ani_long_t_t[ani_long_t_t < 95] <- NA
display_numbers_matrix <- ani_long_t_t
display_numbers_matrix[is.na(display_numbers_matrix)] <- ""

pdf("ani_matrix_typestrain_SGB.pdf", width = 8, height = 6)
pheatmap(ani_long_t, 
         cluster_rows = TRUE, 
         cluster_cols = F, 
         fontsize_row = 8, 
         fontsize_col = 8,
         color = heatmap_colors,
         main = "ANI Matrix",
         # 显示ani_long_t_t
         display_numbers = display_numbers_matrix,
         # number_format = "%.3f",
         number_color = "black",
         cellwidth = 20, 
         cellheight = 20,
         # clustering_distance_cols = "manhattan",  # 列聚类的距离计算方法
         # clustering_distance_rows = "manhattan",  # 列聚类的距离计算方法
         border_color = NA,
         treeheight_row = 30,  # 行树状图高度
         treeheight_col = 30   # 列树状图高度
)
dev.off()
# library(ggplot2)
# library(reshape2)
# ani_long_melt <- melt(as.matrix(ani_long_t))
# # ani_long_melt$Var1 <- gsub(".fna", "", ani_long_melt$Var1)
# # ani_long_melt$Var2 <- gsub(".fna", "", ani_long_melt$Var2)
# ggplot(ani_long_melt, aes(Var2, Var1, fill = value)) +
#   geom_tile() +
#   #好看的颜色
#   # scale_fill_gradient(low = "white", high = "red") +
#   #范围是从最大值到最小值
#   scale_fill_gradientn(colours = c("white", "#CE1B1E"), limits = c(80, 100)) +
#   labs(x = "SGB", y = "Genome", title = "ANI Matrix") +
#   theme_minimal() +
#   # 超过95的进行标注
#   geom_text(aes(label = ifelse(value > 95, round(value, 1), "")), color = "black", size = 3) +
#   
#   coord_fixed() +
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(
#     text = element_text(size = 12),
#     axis.text.x = element_text(size = 10,angle = 90, hjust = 1,vjust=0.5), 
#     # axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 10,face = "italic"),
#     # strip.text = element_text(size = 12, face = "bold"),
#     legend.position = "right",
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   ) 
# 
# ggsave(filename = "ani_matrix_typestrain_SGB.pdf", width = 8, height = 6)
# 



# sp vs sp



ani_long_t <- ani_long[grep("F.", row.names(ani_long)), grep("F.", colnames(ani_long))] 
rownames(ani_long_t) <- sapply(strsplit(rownames(ani_long_t), "_"), function(x) paste(x[1], collapse = "_"))
colnames(ani_long_t) <- sapply(strsplit(colnames(ani_long_t), "_"), function(x) paste(x[1], collapse = "_"))
library(pheatmap)
library(RColorBrewer)
heatmap_colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
#只显示大于95的值
ani_long_t_t <- ani_long_t
ani_long_t_t <- apply(ani_long_t_t, 2, function(x) round(x, 1))
ani_long_t_t[ani_long_t_t < 95] <- NA
display_numbers_matrix <- ani_long_t_t
display_numbers_matrix[is.na(display_numbers_matrix)] <- ""
pdf("ani_matrix_typestrain_typestrain.pdf", width = 8, height = 6)
pheatmap(ani_long_t, 
         cluster_rows = TRUE, 
         cluster_cols = T, 
         fontsize_row = 8, 
         fontsize_col = 8,
         color = heatmap_colors,
         main = "ANI Matrix between Faecalibacterium species",
         # 显示ani_long_t_t
         display_numbers = display_numbers_matrix,
         # number_format = "%.3f",
         number_color = "black",
         cellwidth = 20, 
         cellheight = 20,
         border_color = NA,
         treeheight_row = 30,  # 行树状图高度
         treeheight_col = 30   # 列树状图高度
)
dev.off()
# # SGB vs SGB
# 
# ani_long_t <- ani_long[grep("SGB", row.names(ani_long)), grep("SGB", colnames(ani_long))] 
# rownames(ani_long_t) <- sapply(strsplit(rownames(ani_long_t), "_"), function(x) paste(x[1], collapse = "_"))
# colnames(ani_long_t) <- sapply(strsplit(colnames(ani_long_t), "_"), function(x) paste(x[1], collapse = "_"))
# ani_long_t[is.na(ani_long_t)] <- 0
# library(pheatmap)
# library(RColorBrewer)
# heatmap_colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
# #只显示大于95的值
# ani_long_t_t <- ani_long_t
# ani_long_t_t <- apply(ani_long_t_t, 2, function(x) round(x, 1))
# ani_long_t_t[ani_long_t_t < 95] <- NA
# display_numbers_matrix <- ani_long_t_t
# display_numbers_matrix[is.na(display_numbers_matrix)] <- ""
# 
# pheatmap(ani_long_t, 
#          cluster_rows = TRUE, 
#          cluster_cols = T, 
#          fontsize_row = 8, 
#          fontsize_col = 8,
#          color = heatmap_colors,
#          main = "ANI Matrix between Faecalibacterium species",
#          # 显示ani_long_t_t
#          display_numbers = display_numbers_matrix,
#          # number_format = "%.3f",
#          number_color = "black",
#          cellwidth = 20, 
#          cellheight = 20,
#          border_color = NA)





















# 多个队列*** 各个菌种的比较 ----
my_data0 <- get(load('E:/OneDrive - hust.edu.cn/项目/chenlab/01_2023_PPI_cross_validation/01_data/feat_meta/PPI_unknown_cohort_name_WGS_SGB.RData'))
feat_list <- my_data0$feat_list
feat_list <- my_pair.table(feat_list)
Faecalibacterium <- T
Faecalibacterium_s <- T
if(Faecalibacterium==T){
  feat_list_Faecalibacterium <- lapply(feat_list, function(x) {
    if(length(grep('Faecalibacterium', colnames(x)))==1){
      col_t <- colnames(x)[grep('Faecalibacterium', colnames(x))]
      x <- data.frame(x[,grep('Faecalibacterium', colnames(x))], check.names = FALSE,row.names = rownames(x))
      colnames(x)= col_t
    }else{
      x <- x[,grep('Faecalibacterium', colnames(x))]
    }
    return(x)
  })
  
  # Faecalibacterium_c_list
  Faecalibacterium_c_list <- lapply(feat_list_Faecalibacterium, function(x) {
    Faecalibacterium_c <- colnames(x)
    return(Faecalibacterium_c)
  })
  
  Faecalibacterium_c <- unlist(Faecalibacterium_c_list)
  Faecalibacterium_c_u <- unique(Faecalibacterium_c)
  which( Faecalibacterium_c%in% Faecalibacterium_c_u[10])
  Faecalibacterium_c_u
  save(Faecalibacterium_c_u, file = "01_data/Faecalibacterium_c_u.RData")
}else{
  if(Faecalibacterium_s==T){
    
    feat_list_Faecalibacterium <- lapply(feat_list, function(x) {
      colnames(x) <- sapply(strsplit(colnames(x),'\\|'),'[',7)
      if(length(unique(colnames(x)))<ncol(x)){
        message('dumplicated species!')
        species_unique <- unique(colnames(x))
        feat_unique <- as.data.frame(matrix(data=0,nrow = nrow(x),ncol = length(species_unique)))
        colnames(feat_unique) <- species_unique
        rownames(feat_unique) <- rownames(x)
        for (species_t in species_unique) {
          match_t <- colnames(x)==species_t
          if(sum(match_t)>1){
            feat_t <- x[,colnames(x)==species_t]
            feat_unique[,species_t] <- rowSums(feat_t)
          }else{
            feat_unique[,species_t] <- x[,species_t]
          }
        }
        x <- feat_unique
        return(x)
      }
    })
    
  }else{
    feat_list_Faecalibacterium <- feat_list
  }
}

feat_list_Faecalibacterium$CF_PRJNA510445 <- NULL
feat_list_Faecalibacterium <- my_pair.table(feat_list_Faecalibacterium)

meta_list <- lapply(my_data0$meta_list, function(x){x[,c('Group','host_age','sex','BMI','disease_stage','country')]})
print(lapply(meta_list, function(x){table(x[,c('Group')])}))
meta_list_Faecalibacterium <- meta_list[names(feat_list_Faecalibacterium)]
# grep('Faecalibacterium',colnames(feat_t))

# # 计算每个队列的Faecalibacterium_c_u的平均值
# Faecalibacterium_c_mean <- lapply(feat_list_Faecalibacterium, function(x) {
#   x_mean <- colMeans(x, na.rm = TRUE)
#   return(x_mean)
# })
# Faecalibacterium_c_mean_df <- do.call(rbind, Faecalibacterium_c_mean)
# 
# # 计算每个队列 meta 不同 Group  的Faecalibacterium_c_u的平均值
# Faecalibacterium_c_mean_meta <- lapply(names(meta_list_Faecalibacterium)[1:10], function(x) {
#   print(paste0("Processing: ", x))
#   meta <- meta_list_Faecalibacterium[[x]]
#   feat <- feat_list_Faecalibacterium[[x]]
#   meta$Faecalibacterium_c_mean <- rowMeans(feat, na.rm = TRUE)
#   meta_grouped <- aggregate(Faecalibacterium_c_mean ~ Group, data = meta, FUN = mean)
#   return(meta_grouped)
# })
# Faecalibacterium_c_mean_meta_df <- do.call(rbind, Faecalibacterium_c_mean_meta)
# # 绘制Faecalibacterium_c_mean_meta_df的柱状图
i=1
feat_list_Faecalibacterium_cohort <- lapply(feat_list_Faecalibacterium, function(x) {
  x$cohort <- names(feat_list_Faecalibacterium)[i]
  i <<- i + 1
  return(x)
})
feat_list_Faecalibacterium_df <- do.call(rbind, feat_list_Faecalibacterium_cohort)
meta_list_Faecalibacterium_df <- do.call(rbind, meta_list_Faecalibacterium)
intersect_rownames <- intersect(rownames(feat_list_Faecalibacterium_df), rownames(meta_list_Faecalibacterium_df))
feat_list_Faecalibacterium_df <- cbind(feat_list_Faecalibacterium_df[intersect_rownames,], meta_list_Faecalibacterium_df[intersect_rownames,])

if(Faecalibacterium!=T){
  if(Faecalibacterium_s==T){
    save(feat_list_Faecalibacterium_df,file='D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_all_s.RData')
  }else{
    save(feat_list_Faecalibacterium_df,file='D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_all.RData')
  }
}


# 宽->长
library(dplyr)
library(tidyr)

feat_list_Faecalibacterium_df_long <- feat_list_Faecalibacterium_df[,c(Faecalibacterium_c_u,'cohort','Group')] %>%
  gather(key = "Species", value = "value", -c(cohort,Group))
feat_list_Faecalibacterium_df_long$value <- feat_list_Faecalibacterium_df_long$value * 100
feat_list_Faecalibacterium_df_long <- feat_list_Faecalibacterium_df_long[grep('s__Faecalibacterium_prausnitzii',feat_list_Faecalibacterium_df_long$Species),]
# unique(Faecalibacterium_meta_df$Species)
feat_list_Faecalibacterium_df_long$Species <- gsub("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|", "", 
                                                   feat_list_Faecalibacterium_df_long$Species,fixed = TRUE)
unique(feat_list_Faecalibacterium_df_long$Species)
feat_list_Faecalibacterium_df_long$Species[grep('15316',feat_list_Faecalibacterium_df_long$Species)] <- 'F.longum'
feat_list_Faecalibacterium_df_long$Species[grep('15332',feat_list_Faecalibacterium_df_long$Species)] <- 'F.prausnitzii'
feat_list_Faecalibacterium_df_long$Species[grep('15318',feat_list_Faecalibacterium_df_long$Species)] <- 'F.duncaniae'
feat_list_Faecalibacterium_df_long$Species[grep('15342',feat_list_Faecalibacterium_df_long$Species)] <- 'F.faecis_F.taiwanense'
feat_list_Faecalibacterium_df_long$Species[grep('15322',feat_list_Faecalibacterium_df_long$Species)] <- 'F.hattorii'
unique(feat_list_Faecalibacterium_df_long$Species)
save(feat_list_Faecalibacterium_df_long, file = "01_data/feat_list_Faecalibacterium_df_long.RData")

# 计算每个队列 不同 Group  不同 cohort的species的平均值
Faecalibacterium_meta_df <- feat_list_Faecalibacterium_df_long %>%
  group_by(cohort, Group, Species) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup()
# Faecalibacterium_meta_df <- Faecalibacterium_meta_df[grep('s__Faecalibacterium_prausnitzii',Faecalibacterium_meta_df$Species),]
# # unique(Faecalibacterium_meta_df$Species)
# Faecalibacterium_meta_df$Species <- gsub("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|", "", 
#                                          Faecalibacterium_meta_df$Species,fixed = TRUE)
# unique(Faecalibacterium_meta_df$Species)
# Faecalibacterium_meta_df$Species[grep('15316',Faecalibacterium_meta_df$Species)] <- 'F.longum'
# Faecalibacterium_meta_df$Species[grep('15332',Faecalibacterium_meta_df$Species)] <- 'F.prausnitzii'
# Faecalibacterium_meta_df$Species[grep('15318',Faecalibacterium_meta_df$Species)] <- 'F.duncaniae'
# Faecalibacterium_meta_df$Species[grep('15342',Faecalibacterium_meta_df$Species)] <- 'F.faecis_F.taiwanense'
# Faecalibacterium_meta_df$Species[grep('15322',Faecalibacterium_meta_df$Species)] <- 'F.hattorii'

# 1) 创建堆积柱状图 ----
feat_list_Faecalibacterium_df_long <- get(load("01_data/feat_list_Faecalibacterium_df_long.RData"))

# 计算每个队列 不同 Group  不同 cohort的species的平均值
Faecalibacterium_meta_df <- feat_list_Faecalibacterium_df_long %>%
  group_by(cohort, Group, Species) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup()

custom_colors <- c("#7FC97F","#BEAED4" ,"#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")

Faecalibacterium_meta_df_h <- Faecalibacterium_meta_df %>%
  filter(Group == 'Control')
# 按照中位数大小排序
feat_list_Faecalibacterium_df_long_h <- feat_list_Faecalibacterium_df_long %>% 
  filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(value)

feat_list_Faecalibacterium_df_long_h <- arrange(feat_list_Faecalibacterium_df_long_h, desc(value))
feat_list_Faecalibacterium_df_long_t <- feat_list_Faecalibacterium_df_long
feat_list_Faecalibacterium_df_long_t$Species <- factor(feat_list_Faecalibacterium_df_long_t$Species, 
                                                       # levels = unique(feat_list_Faecalibacterium_df_long$Species[order(tapply(feat_list_Faecalibacterium_df_long$value, feat_list_Faecalibacterium_df_long$Species, median))])
                                                       levels = feat_list_Faecalibacterium_df_long_h$Species
)
Faecalibacterium_meta_df_h$Species <- factor(Faecalibacterium_meta_df_h$Species, 
                                             levels = feat_list_Faecalibacterium_df_long_h$Species)


Faecalibacterium_meta_df_h_c <- Faecalibacterium_meta_df_h %>% 
  filter(Species == 'F.longum') %>% 
  arrange(desc(mean_value)) %>% 
  select(cohort) %>% 
  ungroup()
Faecalibacterium_meta_df_h_c <- Faecalibacterium_meta_df_h_c$cohort

# Faecalibacterium_meta_df_h_c <- as.character(Faecalibacterium_meta_df_h_c)
Faecalibacterium_meta_df_h$cohort <- factor(Faecalibacterium_meta_df_h$cohort, 
                                            levels = Faecalibacterium_meta_df_h_c)
stacked_barplot <- ggplot(Faecalibacterium_meta_df_h, aes(x = cohort, y =mean_value , fill = Species)) +
  geom_bar(stat = "identity", width = 0.25) +
  # facet_wrap(~ Bacteria_Full, nrow = 1) +
  # facet_grid( ~ gene) +
  theme_minimal() +
  labs(x = "", y = "Abundance %", fill = "Bacteria_Full") +
  theme(
    text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 7),
    axis.text.y = element_text(size = 10),
    # strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # 标注每组的数值
  # geom_text(
  #   aes(label = Count, group = Category),
  #   position = position_stack(vjust = 0.5),
  #   size = 4,
  #   color = "white"
  # ) +
  # # 设置自定义颜色方案
  scale_fill_manual(values = custom_colors)

# 显示图形
print(stacked_barplot)
ggsave(filename = "Faecalibacterium_stacked_bar.pdf",  width = 15, height = 5)

# 1.1) 每个队列的最大值的菌 ----
result <- Faecalibacterium_meta_df_h %>%
  group_by(cohort) %>%          # 按 cohort 分组
  slice_max(mean_value, n = 1)  # 每组内选 mean_value 最大的 1 行
result$Species <- as.character(result$Species)
species_counts <- as.data.frame(table(result$Species))
species_counts$Percentage <- round(species_counts$Freq/sum(species_counts$Freq),2)*100
# pie(species_counts, main = "Species Proportion of Maximum mean_value by Cohort")
colnames(species_counts) <- c('Species','Count','Percentage')
custom_colors <- c("#7FC97F","#BEAED4" ,"#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")
species_counts$Species <- factor(species_counts$Species,levels = c('F.longum','F.faecis_F.taiwanense', 'F.duncaniae', 'F.prausnitzii'))
pie_charts <- ggplot(species_counts, aes(x = "", y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", width = 2) +
  coord_polar("y", start = 0) +
  # facet_grid(Bacteria ~ gene) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(fill = "Species") +
  # 显示具体数目和百分比# paste0(Count, "\n(", round(Percentage, 1), "%)"))
  geom_text(aes(label = paste0(Count, "\n(", round(Percentage, 1), "%)")), 
            position = position_stack(vjust = 0.5),
            size = 4,
            color = "black")

# 显示图形
print(pie_charts) 
ggsave('pie max cohort number.pdf',width = 4,height = 5)

# 2) 丰度比较 ----
feat_list_Faecalibacterium_df_long_h <- feat_list_Faecalibacterium_df_long %>% 
  filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(value)
feat_list_Faecalibacterium_df_long_h <- arrange(feat_list_Faecalibacterium_df_long_h, desc(value))
feat_list_Faecalibacterium_df_long_t <- feat_list_Faecalibacterium_df_long
feat_list_Faecalibacterium_df_long_t$Species <- factor(feat_list_Faecalibacterium_df_long_t$Species, 
                                                       # levels = unique(feat_list_Faecalibacterium_df_long$Species[order(tapply(feat_list_Faecalibacterium_df_long$value, feat_list_Faecalibacterium_df_long$Species, median))])
                                                       levels = feat_list_Faecalibacterium_df_long_h$Species
)

feat_list_Faecalibacterium_df_long_t$Group <- factor(feat_list_Faecalibacterium_df_long_t$Group, levels = c('Control','Case'))
colors <- c('#4A72A6','#E41A1C')

# p1 <- ggboxplot(feat_list_Faecalibacterium_df_long_t, x = "Species", color="black",y = "value",fill = 'Group', 
#                 alpha=0.8,
#                 #na
#                 
#                 palette =c('#4A72A6','#E41A1C'),
#                 # add = "jitter",
#                 width = 0.5,
#                 # outlier.size = 0.5,
#                 outlier.colour = NA
#                 )+ 
#   # geom_hline(yintercept =0.5,color='#dbdcdc')+
#   # geom_hline(yintercept =0.6,color='#ffd09a')+
#   # geom_hline(yintercept =0.7,color='#ffcbd8')+
#   stat_compare_means(aes(group=Group),label = "p.signif",label.y = 10,paired = F)+
#   xlab("16S_genus") + ylab("External AUC")+
#   labs(fill ="Disease type")+
#   # ylim(0.30,0.87)+
#   theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none") 
# p1
p3 <- ggboxplot(feat_list_Faecalibacterium_df_long_t , x = 'Species', y ='value' ,
                fill = 'Group',
                # palette = pal_igv("default")(51)[-c(1,2)],
                width = 0.4,
                # color =  "#386CB0",
                outlier.size = 0.5,
                
                outlier.colour = NA
)+
  
  # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
  scale_fill_manual(values =colors)+
  
  # scale_color_igv(alpha = 0.7)+
  # scale_fill_igv(alpha = 0.3)+
  # scale_y_log10()+
  # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
  #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
  
  stat_compare_means(aes(group=Group),label = "p.signif",label.y = 4,paired = F)+
  # method.args = list(alternative = "less")
  # stat_compare_means(
  #   ref.group = "F.longum",
  #                    label = 'p.signif',
  #                    method='wilcox.test',
  #                    label.y = 4.7,
  #                    paired = T
  #                    # method.args = list(alternative = "less")
  # )+
  # ylim(0,0.06)+
  
  coord_cartesian(ylim = c(-0.001, 0.045)*100)+
  
  ylab("Abundance %")+xlab('')+
  
  # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
  theme_bw() +
  theme(
    legend.position="right",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(angle=90, hjust=1,face = 'italic',size=10),
    axis.text.y=element_text(face = 'plain',size=10),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )


# p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p3)

ggsave(filename = "Faecalibacterium_group.pdf", width = 7, height = 5)

Faecalibacterium_meta_df$Group <- factor(Faecalibacterium_meta_df$Group, levels = c('Control','Case'))
p3 <- ggboxplot(Faecalibacterium_meta_df , x = 'Species', y ='mean_value' ,
                fill = 'Group',
                # palette = pal_igv("default")(51)[-c(1,2)],
                width = 0.4,
                # color =  "#386CB0",
                outlier.size = 0.5
                
                # outlier.colour = NA
)+
  
  # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
  scale_fill_manual(values =colors)+
  
  # scale_color_igv(alpha = 0.7)+
  # scale_fill_igv(alpha = 0.3)+
  # scale_y_log10()+
  # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
  #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
  # stat_compare_means(
  #   ref.group = "F.longum",
  #                    label = 'p.signif',
  #                    method='wilcox.test',
  #                    label.y = 4.7,
  #                    paired = T
  #                    # method.args = list(alternative = "less")
  # )+
  # ylim(0,0.06)+
  coord_cartesian(ylim = c(-0.0001, 0.001)*100)+
  
  ylab("Abundance %")+xlab('')+
  
  # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
  theme_bw() +
  theme(
    legend.position="right",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(angle=90, hjust=1,face = 'italic',size=10),
    axis.text.y=element_text(face = 'plain',size=10),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )


# p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p3)
ggsave(filename = "Faecalibacterium_group_cohort.pdf", width = 5, height = 5)


# no group
# 按照中位数大小排序
feat_list_Faecalibacterium_df_long_h <- feat_list_Faecalibacterium_df_long %>% 
  filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(value)
feat_list_Faecalibacterium_df_long_h <- arrange(feat_list_Faecalibacterium_df_long_h, desc(value))
feat_list_Faecalibacterium_df_long_t <- feat_list_Faecalibacterium_df_long
feat_list_Faecalibacterium_df_long_t$Species <- factor(feat_list_Faecalibacterium_df_long_t$Species, 
                                                       # levels = unique(feat_list_Faecalibacterium_df_long$Species[order(tapply(feat_list_Faecalibacterium_df_long$value, feat_list_Faecalibacterium_df_long$Species, median))])
                                                       levels = feat_list_Faecalibacterium_df_long_h$Species
)

p3 <- ggboxplot(feat_list_Faecalibacterium_df_long_t , x = 'Species', y ='value' ,
                # fill = 'Group',
                # fill = 'Group',
                # palette = pal_igv("default")(51)[-c(1,2)],
                width = 0.4,
                color =  "#386CB0",
                outlier.size = 0.5,
                
                outlier.colour = NA
)+
  
  # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
  scale_fill_manual(values =colors)+
  
  # scale_color_igv(alpha = 0.7)+
  # scale_fill_igv(alpha = 0.3)+
  # scale_y_log10()+
  # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
  #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
  stat_compare_means(
    ref.group = "F.longum",
    label = 'p.signif',
    method='wilcox.test',
    label.y = 3.8,
    paired = T
    # method.args = list(alternative = "less")
  )+
  # ylim(0,0.06)+
  coord_cartesian(ylim = c(-0.001, 0.040)*100)+
  
  ylab("Abundance %")+xlab('')+
  
  # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
  theme_bw() +
  theme(
    legend.position="right",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(angle=90, hjust=1,face = 'italic',size=10),
    axis.text.y=element_text(face = 'plain',size=10),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )


# p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p3)
ggsave(filename = "F.longum.box.pdf", width = 4.5, height = 5)
# ggsave(filename = "F.box.non.pdf", width = 5, height = 7)

# 2.1) prevalence abundance ----

load("E:/OneDrive - hust.edu.cn/项目/chenlab/2025-FP-FD/01_data/feat_list_Faecalibacterium_df_long.RData")
feat_list_Faecalibacterium_df_long <- subset(feat_list_Faecalibacterium_df_long,Group=='Control')
feat_long <- feat_list_Faecalibacterium_df_long[,c("Species","value" )]
feat_long$value <- as.numeric(feat_long$value)

feat_long_h <- feat_long %>% 
  # filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(value)
feat_long_h <- arrange(feat_long_h, desc(value))
feat_long_h_t <- feat_long
feat_long_h_t$Species <- factor(feat_long_h_t$Species, 
                                # levels = unique(feat_list_Faecalibacterium_df_long$Species[order(tapply(feat_list_Faecalibacterium_df_long$value, feat_list_Faecalibacterium_df_long$Species, median))])
                                levels = feat_long_h$Species
)

p3 <- ggboxplot(feat_long_h_t , x = 'Species', y ='value' ,
                # fill = 'Group',
                # fill = 'Group',
                # palette = pal_igv("default")(51)[-c(1,2)],
                width = 0.4,
                color =  "#386CB0",
                outlier.size = 0.5,
                
                outlier.colour = NA
)+
  
  # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
  scale_fill_manual(values =colors)+
  
  # scale_color_igv(alpha = 0.7)+
  # scale_fill_igv(alpha = 0.3)+
  # scale_y_log10()+
  # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
  #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
  stat_compare_means(
    ref.group = "F.longum",
    label = 'p.signif',
    method='wilcox.test',
    label.y = 3.8,
    paired = T
    # method.args = list(alternative = "less")
  )+
  # ylim(0,0.06)+
  # coord_cartesian(ylim = c(-0.001, 0.040)*100)+
  
  ylab("Abundance %")+xlab('')+
  
  # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
  theme_bw() +
  theme(
    legend.position="right",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(angle=90, hjust=1,face = 'italic',size=10),
    axis.text.y=element_text(face = 'plain',size=10),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )


# p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p3)


ggsave(filename = "F.longum.box.control.pdf", width = 4.5, height = 5)
# ggsave(filename = "F.box.non.pdf", width = 5, height = 7)

library(tidyr)
#abundance
feat_long_a <- feat_long %>% 
  # filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(abundance = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(abundance))
#prevalence
feat_long_p <- feat_long %>% 
  group_by(Species) %>%
  summarise(prevalence = sum(value > 0, na.rm = TRUE) / n()) %>%
  ungroup() %>%
  arrange(desc(prevalence))

feat_long_ap <- merge(feat_long_a, feat_long_p, by = "Species")
feat_long_ap <- arrange(feat_long_ap, desc(abundance))
feat_long_ap$Species <- factor(feat_long_ap$Species, levels = feat_long_ap$Species)
feat_long_ap$Species

pA <- ggplot(feat_long_ap, aes(x = Species,  y = abundance)) + #, fill = bacterium
  # geom_bar(stat = "identity",fill='#B3E2CD') + #reorder(bacterium,-mean_abundance)
  geom_bar(stat = "identity",fill='#B3E2CD') + #reorder(bacterium,-mean_abundance)
  ylab("Abundance (mean) %") +
  xlab("") +
  # geom_hline(yintercept = 0.5,color='red')+
  theme(legend.position = 'none')+
  theme_bw() +
  theme(
    legend.position="none",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.y=element_text(face = 'plain',size=10),
    axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1,face = 'italic',size=10),
    # axis.text.x=element_blank(),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )
pA
pP <- ggplot(feat_long_ap, aes(x = Species, y = prevalence*100)) + #fill = bacterium
  # geom_bar(stat = "identity", fill = "#FDCDAC") + #x = reorder(bacterium,-prevalence)
  geom_bar(stat = "identity", fill = "#C9E3FF") + #x = reorder(bacterium,-prevalence)
  ylab("Prevalence %") +
  xlab("") +
  # geom_hline(yintercept = 25,color='red')+
  # ggtitle("Prevalence of Bacteria") +
  theme(legend.position = 'none')+
  theme_bw() +
  theme(
    legend.position="none",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.y=element_text(face = 'plain',size=10),
    axis.text.x=element_blank(),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )
pP

p <- ggarrange(pP,pA, ncol = 1, nrow = 2,
               # labels = c("Prevalence", "B"),
               common.legend = TRUE, legend = "right",
               heights  = c(1, 1.6))
p
ggsave(filename = "F.longum.prevalence.abundance.control.pdf", plot = p, width = 3.5, height = 5)








# 3) species wilcox 的比较 ----
# wilcox 
load("E:/OneDrive - hust.edu.cn/项目/chenlab/2025-FP-FD/01_data/feat_list_Faecalibacterium_df_long.RData")
feat_list_Faecalibacterium_df_long_t <- feat_list_Faecalibacterium_df_long %>% 
  group_by(cohort,Species) %>%
  summarise(pvalue = wilcox.test(value ~ Group, paired = F)$p.value,
            fold_change = ifelse(mean(value[Group == 'Control'], na.rm = TRUE)>mean(value[Group == 'Case'], na.rm = TRUE) ,
                                 (mean(value[Group == 'Control'], na.rm = TRUE) / mean(value[Group == 'Case'], na.rm = TRUE)),
                                 - (mean(value[Group == 'Case'], na.rm = TRUE) / mean(value[Group == 'Control'], na.rm = TRUE))
                                 
            )
  ) %>%
  ungroup() 

feat_list_Faecalibacterium_df_long_t$fold_change[is.na(feat_list_Faecalibacterium_df_long_t$pvalue)] <- 0
feat_list_Faecalibacterium_df_long_t$pvalue[is.na(feat_list_Faecalibacterium_df_long_t$pvalue)] <- 1
# feat_list_Faecalibacterium_df_long_t$fold_change[is.na(feat_list_Faecalibacterium_df_long_t$fold_change)] <- 0
feat_list_Faecalibacterium_df_long_t$fold_change[is.infinite(feat_list_Faecalibacterium_df_long_t$fold_change)] <- 0

# heatmap 
feat_list_Faecalibacterium_df_long_t1 <- feat_list_Faecalibacterium_df_long_t
feat_list_Faecalibacterium_df_long_t1$pvalue <- NULL
feat_list_Faecalibacterium_df_long_t1 <- 
  feat_list_Faecalibacterium_df_long_t1 %>% 
  spread(key = Species, value = fold_change) %>% 
  as.data.frame()

rownames(feat_list_Faecalibacterium_df_long_t1) <- feat_list_Faecalibacterium_df_long_t1$cohort
feat_list_Faecalibacterium_df_long_t1$cohort <- NULL
library(pheatmap)
feat_list_Faecalibacterium_df_long_t1 <- apply(feat_list_Faecalibacterium_df_long_t1, 2, function(x) {
  a <- as.numeric(x)
  return(x)
})
feat_list_Faecalibacterium_df_long_t1[feat_list_Faecalibacterium_df_long_t1>10] <- 10
# feat_list_Faecalibacterium_df_long_t1[abs(feat_list_Faecalibacterium_df_long_t1)<1] <- 0
feat_list_Faecalibacterium_df_long_t1 <- t(feat_list_Faecalibacterium_df_long_t1)

# heatmap 2
feat_list_Faecalibacterium_df_long_t2 <- feat_list_Faecalibacterium_df_long_t
feat_list_Faecalibacterium_df_long_t2$fold_change <- NULL
feat_list_Faecalibacterium_df_long_t2 <- 
  feat_list_Faecalibacterium_df_long_t2 %>% 
  spread(key = Species, value = pvalue) %>% 
  as.data.frame()

rownames(feat_list_Faecalibacterium_df_long_t2) <- feat_list_Faecalibacterium_df_long_t2$cohort
feat_list_Faecalibacterium_df_long_t2$cohort <- NULL
library(pheatmap)
feat_list_Faecalibacterium_df_long_t2 <- apply(feat_list_Faecalibacterium_df_long_t2, 2, function(x) {
  a <- as.numeric(x)
  return(x)
})
# feat_list_Faecalibacterium_df_long_t2[feat_list_Faecalibacterium_df_long_t2>0.05] <- NA

feat_list_Faecalibacterium_df_long_t2 <- t(feat_list_Faecalibacterium_df_long_t2)
feat_list_Faecalibacterium_df_long_t1[feat_list_Faecalibacterium_df_long_t2>0.05] <- 0
data_mark <- data.frame()
for(i in 1:nrow(feat_list_Faecalibacterium_df_long_t2)){
  for(j in 1:ncol(feat_list_Faecalibacterium_df_long_t2)){
    if(feat_list_Faecalibacterium_df_long_t2[i,j] <= 0.001)
    {
      data_mark[i,j]="***"
    }
    else if(feat_list_Faecalibacterium_df_long_t2[i,j] <= 0.01 && feat_list_Faecalibacterium_df_long_t2[i,j] > 0.001)
    {
      data_mark[i,j]="**"
    }
    else if(feat_list_Faecalibacterium_df_long_t2[i,j] <= 0.05 && feat_list_Faecalibacterium_df_long_t2[i,j] > 0.01)
    {
      data_mark[i,j]="*"
    }
    else
    {
      data_mark[i,j]=""
    }
  }
}
pheatmap(feat_list_Faecalibacterium_df_long_t1, 
         cluster_rows = T,
         cluster_cols = T,
         fontsize_row = 8, 
         fontsize_col = 8,
         breaks = seq(-10, 10, by = 0.2),
         color = colorRampPalette(c("#FC8D62" , "white", "#66C2A5"))(100),
         main = "Fold Change of Faecalibacterium Species",
         cellwidth = 10, 
         cellheight = 10,
         border_color = 'black',
         treeheight_row = 0,  # 行树状图高度
         treeheight_col = 0,   # 列树状图高度
         #增加标注
         display_numbers = data_mark
)
pdf





# 4) species lefse ----
my_data0 <- get(load('E:/OneDrive - hust.edu.cn/项目/chenlab/01_2023_PPI_cross_validation/01_data/feat_meta/PPI_unknown_cohort_name_WGS_SGB.RData'))
meta_list <- lapply(my_data0$meta_list, function(x){x[,c('Group','host_age','sex','BMI','disease_stage','country')]})
print(lapply(meta_list, function(x){table(x[,c('Group')])}))

feat_list <- my_data0$feat_list
feat_list <- my_pair.table(feat_list)
# 合并 dis val
feat_list_dis_val
if(mer_disval==T){
  unique_cohort <- table(paste0(sapply(strsplit(names(feat_list),'_'),'[',1),'_',sapply(strsplit(names(feat_list),'_'),'[',2)))
  disval_cohort <- names(which(unique_cohort==2))
  for (c_t in disval_cohort) {
    feat_list[[c_t]] <- do.call(rbind,feat_list[grep(c_t,names(feat_list))])
    feat_list[grep(paste0(c_t,'_'),names(feat_list))] <- NULL
    meta_list[[c_t]] <- do.call(rbind,meta_list[grep(c_t,names(meta_list))])
    meta_list[grep(paste0(c_t,'_'),names(meta_list))] <- NULL
  }
}


# feat_list_Faecalibacterium <- lapply(feat_list, function(x) {
#   if(length(grep('Faecalibacterium', colnames(x)))==1){
#     col_t <- colnames(x)[grep('Faecalibacterium', colnames(x))]
#     x <- data.frame(x[,grep('Faecalibacterium', colnames(x))], check.names = FALSE,row.names = rownames(x))
#     colnames(x)= col_t
#   }else{
#     x <- x[,grep('Faecalibacterium', colnames(x))]
#   }
#   
#   
#   return(x)
# })
# 
# # Faecalibacterium_c_list
# Faecalibacterium_c_list <- lapply(feat_list_Faecalibacterium, function(x) {
#   Faecalibacterium_c <- colnames(x)
#   return(Faecalibacterium_c)
# })
# 
# Faecalibacterium_c <- unlist(Faecalibacterium_c_list)
# Faecalibacterium_c_u <- unique(Faecalibacterium_c)
# which( Faecalibacterium_c%in% Faecalibacterium_c_u[10])
# Faecalibacterium_c_u
# # feat_list_Faecalibacterium$CF_PRJNA510445 <- NULL
# 
# 
# 
# meta_list_Faecalibacterium <- meta_list[names(feat_list_Faecalibacterium)]
# 
# 
# feat_list <- feat_list_Faecalibacterium
# meta_list <- meta_list_Faecalibacterium
# feat_list <- my_pair.table(feat_list)
# meta_list <- lapply(meta_list, function(x){x[,c('Group','host_age','sex','BMI','disease_stage','country')]})
# print(lapply(meta_list, function(x){table(x[,c('Group')])}))
# 
# feat_list$CD_PRJNA487636 <- NULL
# meta_list$CD_PRJNA487636 <- NULL
# feat_list$HCC_PRJNA932948 <- NULL
# meta_list$HCC_PRJNA932948 <- NULL
# feat_list$ACVD_PRJEB21528 <- NULL
# meta_list$ACVD_PRJEB21528 <- NULL
# feat_list$CRC_PRJNA961076 <- NULL
# meta_list$CRC_PRJNA961076 <- NULL
# # CRC_PRJNA961076
# # MDD_PRJNA762199
# feat_list$MDD_PRJNA762199 <- NULL
# meta_list$MDD_PRJNA762199 <- NULL
# # ASD_PRJEB23052
# feat_list$ASD_PRJEB23052 <- NULL
# meta_list$ASD_PRJEB23052 <- NULL
# PD_PRJNA433459
# feat_list$PD_PRJNA433459 <- NULL
# meta_list$PD_PRJNA433459 <- NULL
# AS_PRJNA375935_dis

marker.adj <- my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,
                            is_plot=F,lda_cutoff=2,nproj_cutoff=1,level='species',change_name=F)
save(marker.adj,file =  paste0('D:/lm_project/FD 和 FP/marker/','marker.SGB.RData'))

load(paste0('D:/lm_project/FD 和 FP/marker/','marker.SGB.RData'))
load("E:/OneDrive - hust.edu.cn/项目/chenlab/2025-FP-FD/01_data/Faecalibacterium_c_u.RData")
marker_data <- marker.adj$marker_data
marker_data <- marker_data[grep('s__Faecalibacterium_prausnitzii',
                                marker_data$scientific_name),]
# adjust origin 
marker_data <- subset(marker_data, 
                      class %in% c("origin"))
marker_data$nrproj <- NULL
marker_data$class <- NULL
marker_data$LDA <- as.numeric(marker_data$LDA)
colnames(marker_data) <- c('Species',"LDA","p_value","cohort"  )
marker_data$Species[grep('15316',marker_data$Species)] <- 'F.longum'
marker_data$Species[grep('15332',marker_data$Species)] <- 'F.prausnitzii'
marker_data$Species[grep('15318',marker_data$Species)] <- 'F.duncaniae'
marker_data$Species[grep('15342',marker_data$Species)] <- 'F.faecis_F.taiwanense'
marker_data$Species[grep('15322',marker_data$Species)] <- 'F.hattorii'
marker_data$Species[grep('t__SGB15317',marker_data$Species)] <- 't__SGB15317'
marker_data$Species[grep('t__SGB15323',marker_data$Species)] <- 't__SGB15323'


feat_list_Faecalibacterium_df_long_t <- marker_data
# heatmap 
feat_list_Faecalibacterium_df_long_t1 <- feat_list_Faecalibacterium_df_long_t
feat_list_Faecalibacterium_df_long_t1$p_value <- NULL
feat_list_Faecalibacterium_df_long_t1 <- 
  feat_list_Faecalibacterium_df_long_t1 %>% 
  spread(key = Species, value = LDA) %>% 
  as.data.frame()
feat_list_Faecalibacterium_df_long_t1[is.na(feat_list_Faecalibacterium_df_long_t1)] <- 0
rownames(feat_list_Faecalibacterium_df_long_t1) <- feat_list_Faecalibacterium_df_long_t1$cohort
feat_list_Faecalibacterium_df_long_t1$cohort <- NULL
library(pheatmap)
feat_list_Faecalibacterium_df_long_t1 <- apply(feat_list_Faecalibacterium_df_long_t1, 2, function(x) {
  a <- as.numeric(x)
  return(x)
})
# feat_list_Faecalibacterium_df_long_t1[feat_list_Faecalibacterium_df_long_t1>10] <- 10
# feat_list_Faecalibacterium_df_long_t1[abs(feat_list_Faecalibacterium_df_long_t1)<1] <- 0
feat_list_Faecalibacterium_df_long_t1 <- t(feat_list_Faecalibacterium_df_long_t1)

# 
# # heatmap 2
# feat_list_Faecalibacterium_df_long_t2 <- feat_list_Faecalibacterium_df_long_t
# feat_list_Faecalibacterium_df_long_t2$fold_change <- NULL
# feat_list_Faecalibacterium_df_long_t2 <- 
#   feat_list_Faecalibacterium_df_long_t2 %>% 
#   spread(key = Species, value = pvalue) %>% 
#   as.data.frame()
# 
# rownames(feat_list_Faecalibacterium_df_long_t2) <- feat_list_Faecalibacterium_df_long_t2$cohort
# feat_list_Faecalibacterium_df_long_t2$cohort <- NULL
# library(pheatmap)
# feat_list_Faecalibacterium_df_long_t2 <- apply(feat_list_Faecalibacterium_df_long_t2, 2, function(x) {
#   a <- as.numeric(x)
#   return(x)
# })
# # feat_list_Faecalibacterium_df_long_t2[feat_list_Faecalibacterium_df_long_t2>0.05] <- NA
# 
# feat_list_Faecalibacterium_df_long_t2 <- t(feat_list_Faecalibacterium_df_long_t2)
# feat_list_Faecalibacterium_df_long_t1[feat_list_Faecalibacterium_df_long_t2>0.05] <- 0
# data_mark <- data.frame()
# for(i in 1:nrow(feat_list_Faecalibacterium_df_long_t2)){
#   for(j in 1:ncol(feat_list_Faecalibacterium_df_long_t2)){
#     if(feat_list_Faecalibacterium_df_long_t2[i,j] <= 0.001)
#     {
#       data_mark[i,j]="***"
#     }
#     else if(feat_list_Faecalibacterium_df_long_t2[i,j] <= 0.01 && feat_list_Faecalibacterium_df_long_t2[i,j] > 0.001)
#     {
#       data_mark[i,j]="**"
#     }
#     else if(feat_list_Faecalibacterium_df_long_t2[i,j] <= 0.05 && feat_list_Faecalibacterium_df_long_t2[i,j] > 0.01)
#     {
#       data_mark[i,j]="*"
#     }
#     else
#     {
#       data_mark[i,j]=""
#     }
#   }
# }
all=T
if(all==T){
  feat_list_Faecalibacterium_df_long_t1 <- cbind(feat_list_Faecalibacterium_df_long_t1, 
                                                 matrix(0, nrow = nrow(feat_list_Faecalibacterium_df_long_t1), 
                                                        ncol = length(
                                                          setdiff(unique(marker.adj$marker_data$project_id), colnames(feat_list_Faecalibacterium_df_long_t1))
                                                        ),
                                                        dimnames = list(rownames(feat_list_Faecalibacterium_df_long_t1),
                                                                        setdiff(unique(marker.adj$marker_data$project_id), colnames(feat_list_Faecalibacterium_df_long_t1))
                                                        )
                                                 )
  )
  #
}
# 首先检查您的R版本和内存使用情况
# sessionInfo()
# memory.limit()
p <- pheatmap(feat_list_Faecalibacterium_df_long_t1, 
              cluster_rows = T,
              cluster_cols = T,
              fontsize_row = 8, 
              fontsize_col = 8,
              breaks = seq(-5, 5, by = 0.2),
              color = colorRampPalette(c( "#66C2A5", "white", "#FC8D62"))(50),
              main = "Lefse analysis of Faecalibacterium Species",
              cellwidth = 10, 
              cellheight = 10,
              border_color = 'black',
              treeheight_row = 10,  # 行树状图高度
              treeheight_col = 10  # 列树状图高度
              #增加标注
              # display_numbers = data_mark
)
pdf('Lefse analysis of Faecalibacterium Species.pdf', width = 8, height = 6)
print(p)
d


# 5) 加上疾病大类 -----
load('E:/OneDrive - hust.edu.cn/项目/chenlab/01_2023_PPI_cross_validation/01_data/coht_disesae.RData')
# load('D:/OneDrive - hust.edu.cn/项目/chenlab/01_2023_PPI_cross_validation/01_data/coht_disesae.RData')
coht_disesae$Disease_type <- as.character(coht_disesae$Disease_type)
coht_disesae$Cohort <- as.character(coht_disesae$Cohort)
coht_disesae$Disease_type[grep('CRC',coht_disesae$Cohort)] <- 'Cancer'
coht_disesae <- rbind(coht_disesae,c('PPI','PPI_self'))
coht_disesae <- rbind(coht_disesae,c('PPI','H2B_self'))
# load('E:/OneDrive - hust.edu.cn/项目/chenlab/01_2023_PPI_cross_validation/01_data/newdata/data_plot.cohort.detail.RData')
# coht_disesae <- unique(data_plot[,c('Disease_type','Cohort')])
# colnames(coht_disesae)[1] <- 'Disease type'
# 
# coht_disesae$`Disease type` <- as.character(coht_disesae$`Disease type`)
# coht_disesae$Cohort <- as.character(coht_disesae$Cohort)
# coht_disesae <- subset(coht_disesae,Cohort%in%intersect(coht_disesae$Cohort,colnames(datas)))
# 
# # coht_disesae <- rbind(coht_disesae,c('PPI','PPI_self'))
# 
# # coht_disesae$`Disease type` <- factor(coht_disesae$`Disease type`,levels = c("PPI","Liver","Cancer","Intestinal","Metabolic","Mental", "Autoimmune"))
# coht_disesae$`Disease type` <- factor(coht_disesae$`Disease type`)
# coht_disesae <- arrange(coht_disesae,`Disease type`,Cohort)
coht_disesae_cohort <- data.frame("cohort_name"=colnames(feat_list_Faecalibacterium_df_long_t1),
                                  "Cohort"=gsub("_dis|_val","",colnames(feat_list_Faecalibacterium_df_long_t1)))
coht_disesae_cohort <- merge(coht_disesae_cohort, coht_disesae, by = "Cohort", all.x = TRUE)
coht_disesae_cohort$Disease_type <- as.character(coht_disesae_cohort$Disease_type)
# colnames(feat_list_Faecalibacterium_df_long_t1)
coht_disesae_cohort$Disease_type <- factor(coht_disesae_cohort$Disease_type,
                                           levels = names(sort(table(coht_disesae_cohort$Disease_type),decreasing = T)))
coht_disesae_cohort <- arrange(coht_disesae_cohort, Disease_type, cohort_name)

row_names <- data.frame(Disease_type=coht_disesae_cohort$Disease_type,row.names = coht_disesae_cohort$cohort_name)

mycolor<-pal_igv("default", alpha = 0.5)(length(as.character(unique(coht_disesae_cohort$Disease_type))))
names(mycolor) <- as.character(unique(coht_disesae_cohort$Disease_type))

ann_colors = list(
  # Disease_type= c(Autoimmune="#8d958f", Intestinal="#BB002199",  Cancer="#cf9198", PPI="#d89c7c", Liver="#aa3e53",  Metabolic="#6e746a",  Mental="#c8ccc1")
  # Disease_type= c(Autoimmune="#8d958f", Intestinal="#BB002199",  Cancer="#cf9198", PPI="#d89c7c", Liver="#aa3e53",  Metabolic="#6e746a",  Mental="#c8ccc1")
  
  Disease_type=mycolor
  
)

p <- ComplexHeatmap::pheatmap(feat_list_Faecalibacterium_df_long_t1[,coht_disesae_cohort$cohort_name], 
                              annotation_col=row_names,
                              annotation_colors = ann_colors,
                              gaps_col = cumsum(table(coht_disesae_cohort$Disease_type)),
                              
                              # gaps_col = cumsum(table(coht_disesae$`Disease type`)),
                              cluster_rows = T,
                              cluster_cols = F,
                              fontsize_row = 8, 
                              fontsize_col = 8,
                              breaks = seq(-5, 5, by = 0.2),
                              color = colorRampPalette(c( "#66C2A5", "white", "#FC8D62"))(50),
                              #color legend
                              heatmap_legend_param = list(
                                title = "LDA"
                              ),
                              main = "Lefse analysis of Faecalibacterium Species",
                              cellwidth = 10, 
                              cellheight = 10,
                              # border_color = 'black',
                              border_color=NA,
                              treeheight_row = 10,  # 行树状图高度
                              # treeheight_col = 10  # 列树状图高度
                              #增加标注
                              # display_numbers = data_mark
                              # 设置字体 Arial
)

pdf('Lefse analysis of Faecalibacterium Species disease type.pdf', width = 8, height = 6)
print(p)
dev.off()

pdf('Lefse analysis of Faecalibacterium Species disease type all.pdf', width = 12, height = 5)
print(p)
dev.off()

pdf('Lefse analysis of Faecalibacterium Species disease type adjust all.pdf', width = 12, height = 5)
print(p)
dev.off()

# 统计正负个数
rowSums(feat_list_Faecalibacterium_df_long_t1>0)
rowSums(feat_list_Faecalibacterium_df_long_t1<0)

# 统计样本数目 疾病队列
coht_disesae_cohort$disease <- sapply(str_split(coht_disesae_cohort$Cohort,'_'),'[',1)
coht_disesae_cohort$project_id <- sapply(str_split(coht_disesae_cohort$Cohort,'_'),'[',2)
library(readxl)
cohort_information <- read_excel("./01_data/cohort information.xlsx",sheet = 1)
cohort_information <- cohort_information[,c("project id","PMID")]
colnames(cohort_information) <- c("project_id","PMID")
cohort_information <- unique(cohort_information)
coht_disesae_cohort_all <- merge(coht_disesae_cohort,cohort_information,all.x = all)
coht_disesae_cohort_all$project_id[coht_disesae_cohort_all$project_id=='self'] <- 'PRJCA016454'
write.csv(coht_disesae_cohort_all,file = 'coht_disesae_cohort_all.csv',quote = F)

# 6) 加上其它菌的比较 ----
{
  load(paste0('D:/lm_project/FD 和 FP/marker/','marker.SGB.RData'))
  load("E:/OneDrive - hust.edu.cn/项目/chenlab/2025-FP-FD/01_data/Faecalibacterium_c_u.RData")
  
  Faecalibacterium=F
  Faecalibacterium_s=T
  if(Faecalibacterium!=T){
    if(Faecalibacterium_s==T){
      feat_list_Faecalibacterium_df <- get(load('D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_all_s.RData'))
    }else{
      feat_list_Faecalibacterium_df <- get(load('D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_all.RData'))
    }
  }
  Faecalibacterium_prausnitzii_s <- data.frame(feat_list_Faecalibacterium_df[,grep('Faecalibacterium_prausnitzii',colnames(feat_list_Faecalibacterium_df))],
                                               row.names = rownames(feat_list_Faecalibacterium_df))
  colnames(Faecalibacterium_prausnitzii_s) <- 
    'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__SGBall'
  
  Faecalibacterium_s=F
  if(Faecalibacterium!=T){
    if(Faecalibacterium_s==T){
      feat_list_Faecalibacterium_df <- get(load('D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_all_s.RData'))
    }else{
      feat_list_Faecalibacterium_df <- get(load('D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_all.RData'))
    }
  }
  Faecalibacterium_prausnitzii_s <- Faecalibacterium_prausnitzii_s[rownames(feat_list_Faecalibacterium_df),]
  feat_list_Faecalibacterium_df <- cbind(feat_list_Faecalibacterium_df,
                                         'k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__SGBall'=
                                           Faecalibacterium_prausnitzii_s)
  
  
  # 宽->长
  library(dplyr)
  library(tidyr)
  feat_list_Faecalibacterium_df$host_age <- NULL
  feat_list_Faecalibacterium_df$sex <- NULL
  feat_list_Faecalibacterium_df$BMI <- NULL
  feat_list_Faecalibacterium_df$disease_stage <- NULL
  feat_list_Faecalibacterium_df$country <- NULL
  
  feat_list_Faecalibacterium_df_long <- feat_list_Faecalibacterium_df %>%
    gather(key = "Species", value = "value", -c(cohort,Group))
  feat_list_Faecalibacterium_df_long$value <- feat_list_Faecalibacterium_df_long$value * 100
  
  # # feat_list_Faecalibacterium_df_long <- feat_list_Faecalibacterium_df_long[grep('s__Faecalibacterium_prausnitzii',feat_list_Faecalibacterium_df_long$Species),]
  # # # unique(Faecalibacterium_meta_df$Species)
  # # feat_list_Faecalibacterium_df_long$Species <- gsub("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|", "", 
  # #                                                    feat_list_Faecalibacterium_df_long$Species,fixed = TRUE)
  # # unique(feat_list_Faecalibacterium_df_long$Species)
  # feat_list_Faecalibacterium_df_long$Species[grep('15316',feat_list_Faecalibacterium_df_long$Species)] <- 'F.longum'
  # feat_list_Faecalibacterium_df_long$Species[grep('15332',feat_list_Faecalibacterium_df_long$Species)] <- 'F.prausnitzii'
  # feat_list_Faecalibacterium_df_long$Species[grep('15318',feat_list_Faecalibacterium_df_long$Species)] <- 'F.duncaniae'
  # feat_list_Faecalibacterium_df_long$Species[grep('15342',feat_list_Faecalibacterium_df_long$Species)] <- 'F.faecis_F.taiwanense'
  # feat_lt_Faecalibacterium_df_long$Species[grep('15322',feat_list_Faecalibacterium_df_long$Species)] <- 'F.hattorii'
  # # unique(feat_list_Faecalibacterium_df_long$Species)
  # save(feat_list_Faecalibacterium_df_long, file = "D:/lm_project/FD 和 FP/01_data/feat_list_Faecalibacterium_df_long_all.RData")
  
  # 计算每个队列 不同 Group  不同 cohort的species的平均值
  Faecalibacterium_meta_df <- feat_list_Faecalibacterium_df_long %>%
    group_by(Group, Species) %>% #cohort
    summarise(mean_value = mean(value, na.rm = TRUE),
              median_value = median(value, na.rm = TRUE)) %>%
    ungroup()
  # Faecalibacterium_meta_df <- Faecalibacterium_meta_df[grep('s__Faecalibacterium_prausnitzii',Faecalibacterium_meta_df$Species),]
  # # unique(Faecalibacterium_meta_df$Species)
  # Faecalibacterium_meta_df$Species <- gsub("k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|", "", 
  #                                          Faecalibacterium_meta_df$Species,fixed = TRUE)
  # unique(Faecalibacterium_meta_df$Species)
  # Faecalibacterium_meta_df$Species[grep('15316',Faecalibacterium_meta_df$Species)] <- 'F.longum'
  # Faecalibacterium_meta_df$Species[grep('15332',Faecalibacterium_meta_df$Species)] <- 'F.prausnitzii'
  # Faecalibacterium_meta_df$Species[grep('15318',Faecalibacterium_meta_df$Species)] <- 'F.duncaniae'
  # Faecalibacterium_meta_df$Species[grep('15342',Faecalibacterium_meta_df$Species)] <- 'F.faecis_F.taiwanense'
  # Faecalibacterium_meta_df$Species[grep('15322',Faecalibacterium_meta_df$Species)] <- 'F.hattorii'
  
  Faecalibacterium_meta_df_Group <- subset(Faecalibacterium_meta_df,Group=='Control')
  Faecalibacterium_meta_df_Group <- arrange(Faecalibacterium_meta_df_Group,desc(mean_value))
  Faecalibacterium_meta_df_Group$mean_value_index <- 1:nrow(Faecalibacterium_meta_df_Group)
  Faecalibacterium_meta_df_Group <- arrange(Faecalibacterium_meta_df_Group,desc(median_value))
  Faecalibacterium_meta_df_Group$median_value_index <- 1:nrow(Faecalibacterium_meta_df_Group)
  # Faecalibacterium_meta_df_Group1 <- Faecalibacterium_meta_df_Group[1:20,]
  Faecalibacterium_meta_df_Group <- rbind(Faecalibacterium_meta_df_Group[2:6,],
                                          Faecalibacterium_meta_df_Group[grep('Faecalibacterium_prausnitzii',Faecalibacterium_meta_df_Group$Species),])
  
  feat_list_Faecalibacterium_df_long <- subset(feat_list_Faecalibacterium_df_long,Species%in% Faecalibacterium_meta_df_Group$Species )
  feat_list_Faecalibacterium_df_long <- merge(feat_list_Faecalibacterium_df_long, Faecalibacterium_meta_df_Group)
  feat_list_Faecalibacterium_df_long$Species <- paste0(sapply(str_split(feat_list_Faecalibacterium_df_long$Species,'\\|'),'[',7),'_',
                                                       sapply(str_split(feat_list_Faecalibacterium_df_long$Species,'\\|'),'[',8))
  feat_list_Faecalibacterium_df_long$Species[grep('15316',feat_list_Faecalibacterium_df_long$Species)] <- 'F.longum'
  feat_list_Faecalibacterium_df_long$Species[grep('15332',feat_list_Faecalibacterium_df_long$Species)] <- 'F.prausnitzii'
  feat_list_Faecalibacterium_df_long$Species[grep('15318',feat_list_Faecalibacterium_df_long$Species)] <- 'F.duncaniae'
  feat_list_Faecalibacterium_df_long$Species[grep('15342',feat_list_Faecalibacterium_df_long$Species)] <- 'F.faecis_F.taiwanense'
  feat_list_Faecalibacterium_df_long$Species[grep('15322',feat_list_Faecalibacterium_df_long$Species)] <- 'F.hattorii'
  feat_list_Faecalibacterium_df_long <- arrange(feat_list_Faecalibacterium_df_long,(feat_list_Faecalibacterium_df_long$median_value_index))
  
  
  # feat_list_Faecalibacterium_df_long$Species <- paste0(sapply(str_split(feat_list_Faecalibacterium_df_long$Species,'_t__'),'[',1))
  feat_list_Faecalibacterium_df_long$Species <- gsub('s__','',feat_list_Faecalibacterium_df_long$Species)
  feat_list_Faecalibacterium_df_long$Species <- gsub('_|t__',' ',feat_list_Faecalibacterium_df_long$Species)
  feat_list_Faecalibacterium_df_long$Species <- factor(feat_list_Faecalibacterium_df_long$Species,levels = rev(unique(feat_list_Faecalibacterium_df_long$Species)))
}

p6 <- ggplot(feat_list_Faecalibacterium_df_long, aes(x = value, y = Species, fill = median_value)) +
  scale_fill_continuous(low = "white", high = "#66A61E")+
  geom_boxplot(
    outlier.colour = NA
    
  ) +
  labs(
    x = "Avundance (100%)", # x轴标签
    y = "Species" # y轴标签
    # title = "RF 10 times 10 folds" # 图表标题
  )+
  theme_bw()+
  xlim(0,30)+
  theme(axis.text.y=element_text( hjust=1,face = 'italic',size=10),
        axis.text.x=element_text(face = 'plain',size=10),
        text = element_text(size=12,face = 'plain',family ='',colour = 'black'),
        legend.position = 'none')
p6

save(feat_list_Faecalibacterium_df_long,file = '01_data/feat_list_Faecalibacterium_df_long.RData')
ggsave(filename = paste0('02_figure/top_20_species_in_health.pdf'),width = 5,height = 4.5)
ggsave(filename = paste0('02_figure/top_species_in_health_top_F.pdf'),width = 5,height = 4.5)



# BGI SNP ----
norm=F
snp_data <- read.table('01_data/abo.fut2.ped', header = F, sep = ' ')
# 检测存在0的行
snp_data <- snp_data[!apply(snp_data[,-c(1:6)], 1, function(x) any(x == 0)), ]
snp_data_map <- read.table('01_data/abo.fut2.map', header = F, sep = '\t')
snp_data_match <- read.table('01_data/snp.rs.txt', header = F, sep = ' ')
abun <- read.table('01_data/Faecalibacterium.mp4.abundance', header = T, sep = '\t')
# abun$k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium <- NULL
# abun$IID <- NULL
# abun$k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii <- NULL
abun <- abun[,-grep('s__Faecalibacterium_SGB',colnames(abun))]
abun_FP <- abun[,-grep('t__|FID|IID', colnames(abun))]
colnames(abun_FP) <- c('g__Faecalibacterium','g__Faecalibacterium.s__Fp')
abun <- abun[,grep('t__|FID', colnames(abun))]
colnames(abun)[2:ncol(abun)] <- sapply(strsplit(colnames(abun)[2:ncol(abun)], "\\."), '[',8)
colnames(abun)[grep('15316',colnames(abun))] <- 'F.longum'
colnames(abun)[grep('15332',colnames(abun))] <- 'F.prausnitzii'
colnames(abun)[grep('15318',colnames(abun))] <- 'F.duncaniae'
colnames(abun)[grep('15342',colnames(abun))] <- 'F.faecis_F.taiwanense'
colnames(abun)[grep('15322',colnames(abun))] <- 'F.hattorii'
abun <- cbind( abun,abun_FP)

print(paste0('微生物丰度样本共',nrow(abun), '个, ',
             'SNP数据样本共', nrow(snp_data), '个, ',
             '取样本集交集共有', length(intersect(snp_data$V2, abun$FID)), '个'))
sample_names_intersect <- intersect(snp_data$V2, abun$FID)
snp_data <- snp_data[snp_data$V2 %in% sample_names_intersect, ]
abun <- abun[abun$FID %in% sample_names_intersect, ]


# 统计组合出现次数和类型
# rs8176719 <- table(paste(snp_data$V7, snp_data$V8, sep = "_")) %>% as.data.frame()
# rs635634 <- table(paste(snp_data$V9, snp_data$V10, sep = "_")) %>% as.data.frame()
# rs7030248 <- table(paste(snp_data$V11, snp_data$V12, sep = "_")) %>% as.data.frame()
# rs1047781 <- table(paste(snp_data$V13, snp_data$V14, sep = "_")) %>% as.data.frame()
# 
# rs8176719
# rs635634
# rs7030248
# rs1047781
# 
# three_group <- table(paste(snp_data$V7, snp_data$V8, 
#                            snp_data$V9, snp_data$V10, 
#                            snp_data$V11, snp_data$V12, sep = "_")) %>% as.data.frame()
# 

blood_df <- data.frame(
  sample = snp_data$V2,
  rs8176719 = paste(snp_data$V7, snp_data$V8, sep = "_"),
  rs635634 = paste(snp_data$V9, snp_data$V10, sep = "_"),
  rs7030248 = paste(snp_data$V11, snp_data$V12, sep = "_"),
  rs1047781 = paste(snp_data$V13, snp_data$V14, sep = "_")
)
blood_df$rs8176719 <- ifelse(blood_df$rs8176719 == "T_T", "0",
                             ifelse(blood_df$rs8176719 == "TC_T", "1",
                                    ifelse(blood_df$rs8176719 == "TC_TC", "2", "err")))
# blood_df$rs635634 <- ifelse(blood_df$rs635634 == "T_T", "0",
#                             ifelse(blood_df$rs635634 == "T_C", "1",
#                                    ifelse(blood_df$rs635634 == "C_C", "2", "err")))
blood_df$rs635634 <- ifelse(blood_df$rs635634 == "C_C", "0",
                            ifelse(blood_df$rs635634 == "T_C", "1",
                                   ifelse(blood_df$rs635634 == "T_T", "2", "err")))
blood_df$rs7030248 <- ifelse(blood_df$rs7030248 == "G_G", "0",
                             ifelse(blood_df$rs7030248 == "A_G", "1",
                                    ifelse(blood_df$rs7030248 == "A_A", "2", "err")))

# blood_df$rs8176719 <- ifelse(blood_df$rs8176719 == "TC_TC", "0", 
#                              ifelse(blood_df$rs8176719 == "TC_T", "1", 
#                                     ifelse(blood_df$rs8176719 == "T_T", "2", "err")))
# blood_df$rs635634 <- ifelse(blood_df$rs635634 == "T_T", "0",
#                                  ifelse(blood_df$rs635634 == "T_C", "1",
#                                         ifelse(blood_df$rs635634 == "C_C", "2", "err")))
# blood_df$rs7030248 <- ifelse(blood_df$rs7030248 == "A_A", "0",
#                                  ifelse(blood_df$rs7030248 == "A_G", "1",
#                                         ifelse(blood_df$rs7030248 == "G_G", "2", "err")))
# 
blood_df$rs1047781 <- ifelse(blood_df$rs1047781 == "A_A", "FUT2 secretors",
                             ifelse(blood_df$rs1047781 == "T_A", "FUT2 secretors",
                                    ifelse(blood_df$rs1047781 == "T_T", "FUT2 non−secretors", "err")))
library(readxl)
blood_type <- read_excel("./01_data/blood.xlsx",sheet = 1)
blood_type_df <- merge(blood_df, blood_type, all.x = TRUE)
rownames(blood_type_df) <- blood_type_df$sample
blood_type_df$Group <- ifelse(blood_type_df$Type%in% c("A","AB"), "A/AB", "B/O")
blood_type_df2 <- blood_type_df
table(blood_type_df$Type)
(table(blood_type_df$Type)/nrow(blood_type_df))*100
table(blood_type_df$rs1047781)

abun$sample <- abun$FID
abun0 <- abun
abun <- abun0
abun$FID <- NULL
abun_t <- abun[,-ncol(abun)]
# abun_t <- abun_t+ 1e-6  # 避免对数变换时出现负无穷大
# abun_t <- apply(abun_t,1,function(x){log(x)-mean(log(x))})
# library(clr)
# abun_t <- apply(abun_t,1,function(x){log(x/geometric.mean(x))})
# abun_t <- apply(abun_t,1,function(x){log(x)})
# 归一化
if(norm==T){
  abun_t <- apply(abun_t, 1, function(x) {
    x <- x / sum(x, na.rm = TRUE) * 100  # 转换为百分比
    return(x)
  })
}else{
  abun_t <- t(abun_t)  # 转置为样本为行，物种为列
}
abun_t[is.na(abun_t) ]<- 0  # 将 NA 替换为 0
abun_t <- t(abun_t)  # 转置回原来的格式
# abun_t <- apply(abun_t, 2, as.numeric)
colnames(abun_t) <- colnames(abun)[-ncol(abun)]
abun <- cbind('sample'=abun[,ncol(abun)], abun_t)


blood_type_abun_df <- merge(blood_type_df, abun, all.x = TRUE)


colors <- c("#1F78B4","#A6CEE3"  )

# species <- 'F.longum'
# species <- 'F.prausnitzii'
# species <- 'F.duncaniae'
# colnames(blood_type_abun_df)[grep(species, colnames(blood_type_abun_df))] <- 'value'



# library(ggpubr,warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
# stat.test <- compare_means(
#   value~Group,data = blood_type_abun_df,
#   method = "wilcox.test"
# )%>%mutate(y.position = seq(max(blood_type_abun_df$value+0.2),max(blood_type_abun_df$value+0.2)*3,length.out=1))
# x=stat.test$p.adj
# stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
# stat.test$p <- round(stat.test$p,3)
# alpha$Group <- factor(alpha$Group,
blood_type_abun_df$rs1047781 <- factor(blood_type_abun_df$rs1047781,levels = c("FUT2 secretors","FUT2 non−secretors"))
# range(blood_type_abun_df$value)
blood_type_abun_df[,-c(1:7)] <- apply(blood_type_abun_df[,-c(1:7)], 2, as.numeric)
#
blood_type_abun_df2 <- blood_type_abun_df
p <- list()
for (species in unique(colnames(blood_type_abun_df)[-c(1:7)])) {
  p[[species]] <- ggboxplot(blood_type_abun_df , x = 'Group', y =species ,
                            fill = 'rs1047781',
                            # palette = pal_igv("default")(51)[-c(1,2)],
                            width = 0.4,
                            # color =  "#386CB0",
                            outlier.size = 0.5,
                            
                            outlier.colour = NA
  )+
    facet_grid(~rs1047781, scales = 'free') +
    # rs1047781
    # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
    # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
    # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
    scale_fill_manual(values =colors)+
    # stat_pvalue_manual(stat.test,label = "p")+
    # scale_color_igv(alpha = 0.7)+
    # scale_fill_igv(alpha = 0.3)+
    # scale_y_log10()+
    # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
    #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
    # stat_compare_means(
    #   ref.group = "F.longum",
    #                    label = 'p.signif',
    #                    method='wilcox.test',
    #                    label.y = 4.7,
    #                    paired = T
    #                    # method.args = list(alternative = "less")
    # )+
    stat_compare_means( 
      label = 'p.signif',
      #                    method='wilcox.test',
      # label.y = 6,
      method='wilcox.test'
      # method.args = list(alternative = "less")
    )+
    # ylim(0,0.06)+
    # coord_cartesian(ylim = c(-0.001, 0.08)*100)+
    
    ylab("Abundance %")+xlab('')+
    
    # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
    theme_bw() +
    theme(
      legend.position="none",
      # axis.text.x=element_blank(),
      # axis.ticks.x = element_blank(),
      axis.text.x=element_text( hjust=1,face = 'plain',size=10),
      axis.text.y=element_text(face = 'plain',size=10),
      text = element_text(size=12,face = 'plain',family ='',colour = 'black')
    )+
    ggtitle(species)
}
P <- ggarrange(p$F.longum, 
               p$F.prausnitzii,
               p$F.faecis_F.taiwanense, 
               p$F.duncaniae, 
               p$F.hattorii,
               p$t__SGB15323,
               p$t__SGB15317,
               p$t__SGB15326,
               p$t__SGB15333,
               ncol = 3,
               nrow = 3)
print(P)


P <- ggarrange(p$F.longum, 
               p$F.prausnitzii,
               p$F.faecis_F.taiwanense, 
               p$F.duncaniae, 
               p$F.hattorii,
               p$t__SGB15323,
               p$t__SGB15317,
               p$t__SGB15326,
               p$t__SGB15333,
               p$g__Faecalibacterium,
               p$g__Faecalibacterium.s__Fp,
               ncol = 4,
               nrow = 3)
print(P)

# library(ggplotify)
# P <- as.ggplot(P)
ggsave(filename = "blood_group_BGI_all.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_nor_norm.pdf", width = 12, height = 12)


# sum的计算
blood_type_abun_df_t <- blood_type_abun_df
blood_type_abun_df_t$genus <- rowSums(blood_type_abun_df_t[,-c(1:7)])  # 计算总丰度
p3 <- ggboxplot(blood_type_abun_df_t , x = 'Group', y ='genus' ,
                fill = 'rs1047781',
                # palette = pal_igv("default")(51)[-c(1,2)],
                width = 0.4,
                # color =  "#386CB0",
                outlier.size = 0.5,
                
                outlier.colour = NA
)+
  facet_grid(~rs1047781, scales = 'free') +
  # rs1047781
  # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
  scale_fill_manual(values =colors)+
  # stat_pvalue_manual(stat.test,label = "p")+
  # scale_color_igv(alpha = 0.7)+
  # scale_fill_igv(alpha = 0.3)+
  # scale_y_log10()+
  # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
  #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
  # stat_compare_means(
  #   ref.group = "F.longum",
  #                    label = 'p.signif',
  #                    method='wilcox.test',
  #                    label.y = 4.7,
  #                    paired = T
  #                    # method.args = list(alternative = "less")
  # )+
  stat_compare_means(
    label = 'p.signif',
    #                    method='wilcox.test',
    # label.y = 6,
    method='wilcox.test'
    # method.args = list(alternative = "less")
  )+
  # ylim(0,0.06)+
  # coord_cartesian(ylim = c(-0.001, 0.08)*100)+
  
  ylab("Abundance %")+xlab('')+
  
  # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
  theme_bw() +
  theme(
    legend.position="right",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text( hjust=1,face = 'plain',size=10),
    axis.text.y=element_text(face = 'plain',size=10),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )+
  ggtitle('Sum of all species')

print(p3)
ggsave(filename = "blood_group_BGI_sum.pdf", width = 6, height = 4)

# p3 <- ggboxplot(blood_type_abun_df , x = 'Group', y ='F.longum' ,
#                 fill = 'rs1047781',
#                 # palette = pal_igv("default")(51)[-c(1,2)],
#                 width = 0.4,
#                 # color =  "#386CB0",
#                 outlier.size = 0.5,
#                 
#                 outlier.colour = NA
# )+
#   facet_grid(~rs1047781, scales = 'free') +
#   # rs1047781
#   # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
#   # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
#   # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
#   scale_fill_manual(values =colors)+
#   # stat_pvalue_manual(stat.test,label = "p")+
#   # scale_color_igv(alpha = 0.7)+
#   # scale_fill_igv(alpha = 0.3)+
#   # scale_y_log10()+
#   # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
#   #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
#   # stat_compare_means(
#   #   ref.group = "F.longum",
#   #                    label = 'p.signif',
#   #                    method='wilcox.test',
#   #                    label.y = 4.7,
#   #                    paired = T
#   #                    # method.args = list(alternative = "less")
#   # )+
#   stat_compare_means( 
#     label = 'p.signif',
#     #                    method='wilcox.test',
#     # label.y = 6,
#     method='wilcox.test'
#     # method.args = list(alternative = "less")
#   )+
#   # ylim(0,0.06)+
#   # coord_cartesian(ylim = c(-0.001, 0.08)*100)+
#   
#   ylab("Abundance %")+xlab('')+
#   
#   # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
#   theme_bw() +
#   theme(
#     legend.position="right",
#     # axis.text.x=element_blank(),
#     # axis.ticks.x = element_blank(),
#     axis.text.x=element_text( hjust=1,face = 'plain',size=10),
#     axis.text.y=element_text(face = 'plain',size=10),
#     text = element_text(size=12,face = 'plain',family ='',colour = 'black')
#   )
# 
# # p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")
# 
# print(p3)


# ggsave(filename = "blood_group5.pdf", width = 7, height = 5)
# ggsave(filename = "blood_group6.pdf", width = 6, height = 5)
ggsave(filename = "blood_group10.pdf", width = 6, height = 5)






# BGI SNP new ----

{
  norm=F
  snp_data <- read.table('01_data/FP.new/abo.fut2.ped', header = F, sep = ' ')
  # 检测存在0的行
  snp_data <- snp_data[!apply(snp_data[,-c(1:6)], 1, function(x) any(x == 0)), ]
  snp_data_map <- read.table('01_data/FP.new/abo.fut2.map', header = F, sep = '\t')
  snp_data_match <- read.table('01_data/FP.new/snp.rs.txt', header = F, sep = ' ')
  abun <- read.table('01_data/FP.new/gut.pheno.Faecalibacterium.ave', header = T, sep = '\t')
  # abun$k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium <- NULL
  # abun$IID <- NULL
  # abun$k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii <- NULL
  abun <- abun[,-grep('s__Faecalibacterium_SGB',colnames(abun))]
  abun_FP <- abun[,-grep('t__|FID|IID', colnames(abun))]
  colnames(abun_FP) <- c('g__Faecalibacterium','g__Faecalibacterium.s__Fp')
  abun <- abun[,grep('t__|FID', colnames(abun))]
  colnames(abun)[2:ncol(abun)] <- sapply(strsplit(colnames(abun)[2:ncol(abun)], "\\."), '[',2)
  colnames(abun)[grep('15316',colnames(abun))] <- 'F.longum'
  colnames(abun)[grep('15332',colnames(abun))] <- 'F.prausnitzii'
  colnames(abun)[grep('15318',colnames(abun))] <- 'F.duncaniae'
  colnames(abun)[grep('15342',colnames(abun))] <- 'F.faecis_F.taiwanense'
  colnames(abun)[grep('15322',colnames(abun))] <- 'F.hattorii'
  abun <- cbind( abun,abun_FP)
  
  print(paste0('微生物丰度样本共',nrow(abun), '个, ',
               'SNP数据样本共', nrow(snp_data), '个, ',
               '取样本集交集共有', length(intersect(snp_data$V2, abun$FID)), '个'))
  sample_names_intersect <- intersect(snp_data$V2, abun$FID)
  snp_data <- snp_data[snp_data$V2 %in% sample_names_intersect, ]
  abun <- abun[abun$FID %in% sample_names_intersect, ]
  
  
  # 统计组合出现次数和类型
  rs8176747 <- table(paste(snp_data$V7, snp_data$V8, sep = "_")) %>% as.data.frame()
  rs8176746 <- table(paste(snp_data$V9, snp_data$V10, sep = "_")) %>% as.data.frame()
  rs8176719  <- table(paste(snp_data$V11, snp_data$V12, sep = "_")) %>% as.data.frame()
  rs635634 <- table(paste(snp_data$V13, snp_data$V14, sep = "_")) %>% as.data.frame()
  rs7030248 <- table(paste(snp_data$V15, snp_data$V16, sep = "_")) %>% as.data.frame()
  rs1047781 <- table(paste(snp_data$V17, snp_data$V18, sep = "_")) %>% as.data.frame()
  # 
  rs8176747
  rs8176746
  #
  rs8176719
  rs635634
  rs7030248
  rs1047781
  # three_group <- table(paste(snp_data$V7, snp_data$V8, 
  #                            snp_data$V9, snp_data$V10, 
  #                            snp_data$V11, snp_data$V12, sep = "_")) %>% as.data.frame()
  # 
  
  blood_df <- data.frame(
    sample = snp_data$V2,
    rs8176747= paste(snp_data$V7, snp_data$V8, sep = "_"),
    rs8176746 = paste(snp_data$V9, snp_data$V10, sep = "_"),
    rs8176719 = paste(snp_data$V11, snp_data$V12, sep = "_"),
    rs635634 = paste(snp_data$V13, snp_data$V14, sep = "_"),
    rs7030248 = paste(snp_data$V15, snp_data$V16, sep = "_"),
    rs1047781 = paste(snp_data$V17, snp_data$V18, sep = "_")
  )
  blood_df_supp <- blood_df
  colnames(blood_df_supp) <- c('t_sample',
                               't_rs8176747',
                               't_rs8176746',
                               't_rs8176719',
                               't_rs635634',
                               't_rs7030248',
                               't_rs1047781')
  blood_df <- cbind(blood_df,blood_df_supp)
  
  blood_df$rs8176747 <- ifelse(blood_df$rs8176747 == "G_G", "0",
                               ifelse(blood_df$rs8176747 == "G_C", "1",
                                      ifelse(blood_df$rs8176747 == "C_C", "2", "err")))
  
  blood_df$rs8176719 <- ifelse(blood_df$rs8176719 == "T_T", "0",
                               ifelse(blood_df$rs8176719 == "TC_T", "1",
                                      ifelse(blood_df$rs8176719 == "TC_TC", "2", "err")))
  
  blood_df$rs8176746 <- ifelse(blood_df$rs8176746 == "G_G", "0",
                               ifelse(blood_df$rs8176746 == "T_G", "1",
                                      ifelse(blood_df$rs8176746 == "T_T", "2", "err")))
  
  blood_df$rs635634 <- ifelse(blood_df$rs635634 == "C_C", "0",
                              ifelse(blood_df$rs635634 == "T_C", "1",
                                     ifelse(blood_df$rs635634 == "T_T", "2", "err")))
  blood_df$rs7030248 <- ifelse(blood_df$rs7030248 == "G_G", "0",
                               ifelse(blood_df$rs7030248 == "A_G", "1",
                                      ifelse(blood_df$rs7030248 == "A_A", "2", "err")))
  blood_df$rs1047781_t <- blood_df$rs1047781
  blood_df$rs1047781 <- ifelse(blood_df$rs1047781 == "A_A", "FUT2 secretors",
                               ifelse(blood_df$rs1047781 == "T_A", "FUT2 secretors",
                                      ifelse(blood_df$rs1047781 == "T_T", "FUT2 non−secretors", "err")))
  blood_df$rs1047781_t <- ifelse(blood_df$rs1047781_t == "A_A", "FUT2 secretors AA",
                                 ifelse(blood_df$rs1047781_t == "T_A", "FUT2 secretors Aa",
                                        ifelse(blood_df$rs1047781_t == "T_T", "FUT2 non−secretors aa", "err")))
  
  library(readxl)
  blood_type1 <- read_excel("./01_data/FP.new/blood_all.xlsx",sheet = 1)
  blood_type2 <- read_excel("./01_data/FP.new/blood_all.xlsx",sheet = 2)
  blood_type_df <- merge(blood_df, blood_type1, all.x = TRUE)
  blood_type_df <- merge(blood_type_df, blood_type2, all.x = TRUE)
  # 血型判断 
  print(paste0('总样本数为：',nrow(blood_type_df),
               ', 血型匹配样本数：',sum(blood_type_df$Type1== blood_type_df$Type2),
               ', 血型匹配相似度：',round(sum(blood_type_df$Type1== blood_type_df$Type2)/nrow(blood_type_df)*100,2),'%'))
  
  
  rownames(blood_type_df) <- blood_type_df$sample
  blood_type_df$Group1 <- ifelse(blood_type_df$Type1%in% c("A","AB"), "A/AB", "B/O")
  blood_type_df$Group2 <- ifelse(blood_type_df$Type2%in% c("A","AB"), "A/AB", "B/O")
  
  table(blood_type_df$Type1)
  (table(blood_type_df$Type1)/nrow(blood_type_df))*100
  table(blood_type_df$rs1047781)
  blood_type_df1 <- blood_type_df
  abun$sample <- abun$FID
  abun0 <- abun
  abun <- abun0
  abun$FID <- NULL
  abun_t <- abun[,-c(ncol(abun),ncol(abun)-1,ncol(abun)-2)]
  # abun_t <- abun_t+ 1e-6  # 避免对数变换时出现负无穷大
  # abun_t <- apply(abun_t,1,function(x){log(x)-mean(log(x))})
  # library(clr)
  # abun_t <- apply(abun_t,1,function(x){log(x/geometric.mean(x))})
  # abun_t <- apply(abun_t,1,function(x){log(x)})
  # 归一化
  if(norm==T){
    abun_t <- apply(abun_t, 1, function(x) {
      x <- x / sum(x, na.rm = TRUE) * 100  # 转换为百分比
      return(x)
    })
  }else{
    abun_t <- t(abun_t)  # 转置为样本为行，物种为列
  }
  abun_t[is.na(abun_t) ]<- 0  # 将 NA 替换为 0
  abun_t <- t(abun_t)  # 转置回原来的格式
  # abun_t <- apply(abun_t, 2, as.numeric)
  colnames(abun_t) <- colnames(abun)[-c(ncol(abun),ncol(abun)-1,ncol(abun)-2)]
  abun_t <- cbind(abun_t,abun_FP)
  if(norm!=T){
    abun_t <- abun_t*100
  }
  abun <- cbind('sample'=abun[,ncol(abun)], abun_t)
  
  
  blood_type_abun_df <- merge(blood_type_df, abun, all.x = TRUE)
  
  
  colors <- c("#1F78B4","#A6CEE3"  )
  # colors <- c( "#33A02C","#B2DF8A")
  
  # species <- 'F.longum'
  # species <- 'F.prausnitzii'
  # species <- 'F.duncaniae'
  # colnames(blood_type_abun_df)[grep(species, colnames(blood_type_abun_df))] <- 'value'
  
  
  
  # library(ggpubr,warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
  # stat.test <- compare_means(
  #   value~Group,data = blood_type_abun_df,
  #   method = "wilcox.test"
  # )%>%mutate(y.position = seq(max(blood_type_abun_df$value+0.2),max(blood_type_abun_df$value+0.2)*3,length.out=1))
  # x=stat.test$p.adj
  # stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
  # stat.test$p <- round(stat.test$p,3)
  # alpha$Group <- factor(alpha$Group,
  blood_type_abun_df$rs1047781 <- factor(blood_type_abun_df$rs1047781,levels = c("FUT2 secretors","FUT2 non−secretors"))
  # range(blood_type_abun_df$value)
  col_t <- which(colnames(blood_type_abun_df)=='F.prausnitzii')-1
  blood_type_abun_df[,-c(1:col_t)] <- apply(blood_type_abun_df[,-c(1:col_t)], 2, as.numeric)
  #
  blood_type_abun_df1 <- blood_type_abun_df
  
  p <- list()
  # blood_type_abun_df$Type <- blood_type_abun_df$Type1
  # blood_type_abun_df$Group <- blood_type_abun_df$Group1
  
  blood_type_abun_df$Type <- blood_type_abun_df$Type2
  blood_type_abun_df$Group <- blood_type_abun_df$Group2
  
  blood_type_abun_df$Group3 <- 'Others'
  blood_type_abun_df$Group3 [(blood_type_abun_df$Group== 'A/AB')&(blood_type_abun_df$rs1047781== 'FUT2 secretors')] <- 'A/AB FUT2 secretors'
  
  sp <- unique(colnames(blood_type_abun_df)[-c(1:col_t)])
  sp <- sp[-grep('Type|Group', sp)]
  
  #正态分布检验
  a <- apply(blood_type_abun_df[,-c(1:col_t,ncol(blood_type_abun_df),ncol(blood_type_abun_df)-1,ncol(blood_type_abun_df)-2)], 2, shapiro.test)
  a <- sapply(a, function(x) x$p.value)
  a
  hist(blood_type_abun_df$F.longum, col=' steelblue ') 
  # 基本上不服从正态分布检验
  
  for (species in sp) {
    p[[species]] <- blood_type_abun_df %>%
      # filter(rs1047781 != 'FUT2 secretors') %>%
      # filter(rs1047781 != 'FUT2 non−secretors') %>%
      # filter(rs1047781_t == 'FUT2 secretors AA') %>%
      # filter(rs1047781_t == 'FUT2 secretors Aa') %>%
      # filter(rs1047781_t == 'FUT2 non−secretors aa') %>%
      ggboxplot( x = 'Group', y =species , # x = 'Type1'
                 # ggboxplot( x = 'Group3', y =species , # x = 'Type1'
                 # fill = '#1F78B4',
                 # fill = 'rs1047781',
                 fill = 'Group',
                 # fill = 'Group3',
                 # palette = pal_igv("default")(51)[-c(1,2)],
                 width = 0.4,
                 # color =  "#386CB0",
                 outlier.size = 0.5,
                 
                 outlier.colour = NA
      )+
      facet_grid(~rs1047781, scales = 'free') +
      # facet_grid(~rs1047781, scales = 'free') +
      # rs1047781
      # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
      # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
      # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
      scale_fill_manual(values =colors)+
      # stat_pvalue_manual(stat.test,label = "p")+
      # scale_color_igv(alpha = 0.7)+
      # scale_fill_igv(alpha = 0.3)+
      # scale_y_log10()+
      # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
      #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
      # stat_compare_means(
      #   ref.group = "F.longum",
      #                    label = 'p.signif',
      #                    method='wilcox.test',
      #                    label.y = 4.7,
      #                    paired = T
      #                    # method.args = list(alternative = "less")
      # )+
      stat_compare_means( 
        # comparisons = list(c("A", "AB"),
        #                    c("A", "B"),
        #                    c("A","O"),
        #                    c("AB", "B"),
        #                    c("AB", "O"),
        #                    c("B", "O")
        #                    ),
        label = 'p',# 'p.signif'
        #                    method='wilcox.test',
        # label.y = 6,
        method='wilcox.test',
        # method='t.test',
        method.args = list(alternative = "less") #  greater
      )+
      # ylim(0,5)+
      # coord_cartesian(ylim = c(-0.001, 0.08)*100)+
      
      ylab("Abundance %")+xlab('')+
      # ylab("Percentage %") +    xlab('')+
      
      # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
      theme_bw() +
      theme(
        legend.position="none",
        # axis.text.x=element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x=element_text( hjust=1,face = 'plain',size=10),
        axis.text.y=element_text(face = 'plain',size=10),
        text = element_text(size=12,face = 'plain',family ='',colour = 'black')
      )+
      ggtitle(species)
  }
  # 
  # P <- ggarrange(p$F.longum, 
  #                p$F.prausnitzii,
  #                p$F.faecis_F.taiwanense, 
  #                p$F.duncaniae, 
  #                p$F.hattorii,
  #                p$t__SGB15323,
  #                p$t__SGB15317,
  #                p$t__SGB15326,
  #                p$t__SGB15333,
  #                # p$g__Faecalibacterium,
  #                # p$g__Faecalibacterium.s__Fp,
  #                ncol = 3,
  #                nrow = 3)
  
  P <- ggarrange(p$F.longum, 
                 p$F.prausnitzii,
                 p$F.faecis_F.taiwanense, 
                 p$F.duncaniae, 
                 p$F.hattorii,
                 p$t__SGB15323,
                 p$t__SGB15317,
                 p$t__SGB15326,
                 p$t__SGB15333,
                 p$g__Faecalibacterium,
                 p$g__Faecalibacterium.s__Fp,
                 ncol = 4,
                 nrow = 3)
  print(P)
  
}



if(norm!=T){
  save(blood_type_abun_df,file = '01_data/blood_type_abun_df.norm.F.RData')
}
if(norm!=T){
  save(blood_type_abun_df,file = '01_data/blood_type_abun_df.norm.F1.RData')
}
write.csv(blood_type_abun_df,file = 'blood_FUT2_1.csv',quote = F)

# library(ggplotify)
# P <- as.ggplot(P)
ggsave(filename = "blood_group_BGI_all_newdata.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_nor_norm_newdata.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_t_test_less.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_wilcox_less.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_wilcox_less_percentage.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_wilcox_percentage.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_wilcox_percentage_less.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_wilcox_percentage_group3_less.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_t_test.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_t_test_less_type1.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_t_test_less_Group3.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_wilcox_Group3.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type_FUT2all.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type_FUT2_non_secretors.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type_FUT2_secretors.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type_FUT2_secretors_AA.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type_FUT2_secretors_A_a.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_BGI_all_newdata_p_genus_all_type_FUT2_non_secretors_aa.pdf", width = 12, height = 12)
# THSBC SNP new ----

{ 
  
  norm=F
  library(openxlsx)
  
  abun_snp = read.xlsx('./01_data/share_hust/share_hust/extract_faecali_and_snps_thsbc.xlsx')
  faecali_codebook = read.xlsx('./01_data/share_hust/share_hust/faecali_codebook.xlsx')
  snp_codebook = read.xlsx('./01_data/share_hust/share_hust/snp_codebook.xlsx')
  # table(abun_snp$batch)
  colnames(abun_snp)[3:13] <- faecali_codebook$name
  # abun_snp[,3:13] <- apply(abun_snp[,3:13],2,as.numeric)
  # colnames(abun_snp)[14:19] <- paste0(sapply(str_split(colnames(abun_snp)[14:19],':|_'),'[',2))
  colnames(abun_snp)[14:19] <- snp_codebook$SNP
  # table(abun_snp$rs1047781)
  # table(abun_snp$rs7030248)
  # table(abun_snp$rs635634)
  # table(abun_snp$rs8176719)
  # table(abun_snp$rs8176746)
  # table(abun_snp$rs8176747)
  blood_df <- abun_snp
  blood_df$rs1047781_t <- blood_df$rs1047781
  blood_df$rs1047781 <- ifelse(blood_df$rs1047781 == 0, "FUT2 secretors",
                               ifelse(blood_df$rs1047781 == 1, "FUT2 secretors",
                                      ifelse(blood_df$rs1047781 == 2, "FUT2 non−secretors", "err")))
  library(readxl)
  blood_type1 <- read_excel("./01_data/FP.new/blood_all.xlsx",sheet = 1)
  blood_type2 <- read_excel("./01_data/FP.new/blood_all.xlsx",sheet = 2)
  blood_type_df <- merge(blood_df, blood_type1, all.x = TRUE)
  blood_type_df <- merge(blood_type_df, blood_type2, all.x = TRUE)
  # 血型判断 
  print(paste0('总样本数为：',nrow(blood_type_df),
               ', 血型匹配样本数：',sum(blood_type_df$Type1== blood_type_df$Type2),
               ', 血型匹配相似度：',round(sum(blood_type_df$Type1== blood_type_df$Type2)/nrow(blood_type_df)*100,2),'%'))
  
  rownames(blood_type_df) <- blood_type_df$id
  blood_type_df$Group1 <- ifelse(blood_type_df$Type1%in% c("A","AB"), "A/AB", "B/O")
  blood_type_df$Group2 <- ifelse(blood_type_df$Type2%in% c("A","AB"), "A/AB", "B/O")
  
  table(blood_type_df$Type1)
  (table(blood_type_df$Type1)/nrow(blood_type_df))*100
  table(blood_type_df$rs1047781)
  
  blood_type_df1 <- blood_type_df
  
  abun <- blood_type_df1[,7:17]
  blood_type_df1[,7:17] <- NULL
  abun_t <- abun[,-c(1,2)]
  abun_FP <- abun[,c(1,2)]
  # abun_t <- abun_t+ 1e-6  # 避免对数变换时出现负无穷大
  # abun_t <- apply(abun_t,1,function(x){log(x)-mean(log(x))})
  # library(clr)
  # abun_t <- apply(abun_t,1,function(x){log(x/geometric.mean(x))})
  # abun_t <- apply(abun_t,1,function(x){log(x)})
  # 归一化
  if(norm==T){
    abun_t <- apply(abun_t, 1, function(x) {
      x <- x / sum(x, na.rm = TRUE) * 100  # 转换为百分比
      return(x)
    })
  }else{
    abun_t <- t(abun_t)  # 转置为样本为行，物种为列
  }
  abun_t[is.na(abun_t) ]<- 0  # 将 NA 替换为 0
  abun_t <- t(abun_t)  # 转置回原来的格式
  # abun_t <- apply(abun_t, 2, as.numeric)
  # colnames(abun_t) <- colnames(abun)[-c(ncol(abun),ncol(abun)-1,ncol(abun)-2)]
  abun_t <- cbind(abun_t,abun_FP)
  # if(norm!=T){
  #   abun_t <- abun_t*100
  # }
  # abun <- cbind('sample'=abun[,ncol(abun)], abun_t)
  blood_type_abun_df <- cbind(blood_type_df1, abun_t)
  blood_type_abun_df$s__Faecalibacterium_SGB15315 <- NULL
  
  colnames(blood_type_abun_df)[grep('g__Faecalibacterium',colnames(blood_type_abun_df))] <- 'g__Faecalibacterium'
  colnames(blood_type_abun_df)[grep('s__Faecalibacterium_prausnitzii',colnames(blood_type_abun_df))] <- 'g__Faecalibacterium.s__Fp'
  colnames(blood_type_abun_df)[grep('15316',colnames(blood_type_abun_df))] <- 'F.longum'
  colnames(blood_type_abun_df)[grep('15332',colnames(blood_type_abun_df))] <- 'F.prausnitzii'
  colnames(blood_type_abun_df)[grep('15318',colnames(blood_type_abun_df))] <- 'F.duncaniae'
  colnames(blood_type_abun_df)[grep('15342',colnames(blood_type_abun_df))] <- 'F.faecis_F.taiwanense'
  colnames(blood_type_abun_df)[grep('15322',colnames(blood_type_abun_df))] <- 'F.hattorii'
  colnames(blood_type_abun_df)[grep('15315',colnames(blood_type_abun_df))] <- 'F.butyricigenerans'
  
  
  colors <- c("#1F78B4","#A6CEE3"  )
  # colors <- c( "#33A02C","#B2DF8A")
  
  # species <- 'F.longum'
  # species <- 'F.prausnitzii'
  # species <- 'F.duncaniae'
  # colnames(blood_type_abun_df)[grep(species, colnames(blood_type_abun_df))] <- 'value'
  
  
  
  # library(ggpubr,warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)
  # stat.test <- compare_means(
  #   value~Group,data = blood_type_abun_df,
  #   method = "wilcox.test"
  # )%>%mutate(y.position = seq(max(blood_type_abun_df$value+0.2),max(blood_type_abun_df$value+0.2)*3,length.out=1))
  # x=stat.test$p.adj
  # stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
  # stat.test$p <- round(stat.test$p,3)
  # alpha$Group <- factor(alpha$Group,
  blood_type_abun_df$rs1047781 <- factor(blood_type_abun_df$rs1047781,levels = c("FUT2 secretors","FUT2 non−secretors"))
  # range(blood_type_abun_df$value)
  col_t <- which(colnames(blood_type_abun_df)=='F.longum')-1
  blood_type_abun_df[,-c(1:col_t)] <- apply(blood_type_abun_df[,-c(1:col_t)], 2, as.numeric)
  #
  # blood_type_abun_df1 <- blood_type_abun_df
  
  
  # blood_type_abun_df$Type <- blood_type_abun_df$Type1
  # blood_type_abun_df$Group <- blood_type_abun_df$Group1
  
  blood_type_abun_df$Type <- blood_type_abun_df$Type2
  blood_type_abun_df$Group <- blood_type_abun_df$Group2
  
  blood_type_abun_df$Group3 <- 'Others'
  blood_type_abun_df$Group3 [(blood_type_abun_df$Group== 'A/AB')&(blood_type_abun_df$rs1047781== 'FUT2 secretors')] <- 'A/AB FUT2 secretors'
  
  sp <- unique(colnames(blood_type_abun_df)[-c(1:col_t)])
  sp <- sp[-grep('Type|Group', sp)]
  
  #正态分布检验
  a <- apply(blood_type_abun_df[,-c(1:col_t,ncol(blood_type_abun_df),ncol(blood_type_abun_df)-1,ncol(blood_type_abun_df)-2)], 2, shapiro.test)
  a <- sapply(a, function(x) x$p.value)
  a
  hist(blood_type_abun_df$F.longum, col=' steelblue ') 
  # 基本上不服从正态分布检验
  blood_type_abun_df$Group <- factor(blood_type_abun_df$Group,levels = c('A/AB','B/O'))
  blood_type_abun_df$Group3 <- factor(blood_type_abun_df$Group3,levels = c('A/AB FUT2 secretors','Others'))
  p <- list()
  for (species in sp) {
    p[[species]] <- blood_type_abun_df %>%
      # filter(rs1047781 != 'FUT2 secretors') %>%
      # filter(rs1047781 != 'FUT2 secretors') %>%
      # filter(rs1047781 != 'FUT2 non−secretors') %>%
      # filter(rs1047781_t == 'FUT2 secretors AA') %>%
      # filter(rs1047781_t == 'FUT2 secretors Aa') %>%
      # filter(rs1047781_t == 'FUT2 non−secretors aa') %>%
      # ggboxplot( x = 'Group', y =species , # x = 'Type1'
      ggboxplot( x = 'Group3', y =species , # x = 'Type1'
                 # fill = '#1F78B4',
                 # fill = 'rs1047781',
                 # fill = 'Group',
                 fill = 'Group3',
                 # palette = pal_igv("default")(51)[-c(1,2)],
                 width = 0.4,
                 # color =  "#386CB0",
                 outlier.size = 0.5,
                 
                 outlier.colour = NA
      )+
      # facet_grid(~rs1047781, scales = 'free') +
      # facet_grid(~rs1047781, scales = 'free') +
      # rs1047781
      # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
      # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
      # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
      scale_fill_manual(values =colors)+
      # stat_pvalue_manual(stat.test,label = "p")+
      # scale_color_igv(alpha = 0.7)+
      # scale_fill_igv(alpha = 0.3)+
      # scale_y_log10()+
      # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
      #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
      # stat_compare_means(
      #   ref.group = "F.longum",
      #                    label = 'p.signif',
      #                    method='wilcox.test',
      #                    label.y = 4.7,
      #                    paired = T
      #                    # method.args = list(alternative = "less")
      # )+
      stat_compare_means( 
        # comparisons = list(c("A", "AB"),
        #                    c("A", "B"),
        #                    c("A","O"),
        #                    c("AB", "B"),
        #                    c("AB", "O"),
        #                    c("B", "O")
        #                    ),
        label = 'p',# 'p.signif'
        #                    method='wilcox.test',
        label.y = 6,
        # method='wilcox.test',
        method='t.test',
        method.args = list(alternative = "greater") #  greater
      )+
      # ylim(0,5)+
      # coord_cartesian(ylim = c(-0.001, 0.08)*100)+
      
      ylab("Abundance %")+xlab('')+
      # ylab("Percentage %") +    xlab('')+
      
      # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
      theme_bw() +
      theme(
        legend.position="none",
        # axis.text.x=element_blank(),
        # axis.ticks.x = element_blank(),
        axis.text.x=element_text( hjust=1,face = 'plain',size=10),
        axis.text.y=element_text(face = 'plain',size=10),
        text = element_text(size=12,face = 'plain',family ='',colour = 'black')
      )+
      ggtitle(species)
  }
  # 
  # P <- ggarrange(p$F.longum, 
  #                p$F.prausnitzii,
  #                p$F.faecis_F.taiwanense, 
  #                p$F.duncaniae, 
  #                p$F.hattorii,
  #                p$t__SGB15323,
  #                p$t__SGB15317,
  #                p$t__SGB15326,
  #                p$t__SGB15333,
  #                # p$g__Faecalibacterium,
  #                # p$g__Faecalibacterium.s__Fp,
  #                ncol = 3,
  #                nrow = 3)
  
  P <- ggarrange(p$F.longum, 
                 p$F.prausnitzii,
                 p$F.faecis_F.taiwanense, 
                 p$F.duncaniae, 
                 p$F.hattorii,
                 p$F.butyricigenerans,
                 p$t__SGB15323,
                 p$t__SGB15317,
                 # p$t__SGB15326,
                 # p$t__SGB15333,
                 p$g__Faecalibacterium,
                 p$g__Faecalibacterium.s__Fp,
                 ncol = 4,
                 nrow = 3)
  print(P)
  
}


ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_t_test_greater.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_wilcox_greater.pdf", width = 12, height = 12)

ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_t_test_greater_Group3.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_wilcox_greater_Group3.pdf", width = 12, height = 12)

ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_wilcox_percentage_greater.pdf", width = 12, height = 12)
ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_wilcox_percentage_greater_Group3.pdf", width = 12, height = 12)

table(blood_type_abun_df$Type)
table(blood_type_abun_df$rs1047781)
blood_type_abun_df$rs8176747 <- ifelse(blood_type_abun_df$rs8176747 == "0", "G_G",
                                       ifelse(blood_type_abun_df$rs8176747 == "1", "G_C",
                                              ifelse(blood_type_abun_df$rs8176747 == "2", "C_C", "err")))

blood_type_abun_df$rs8176719 <- ifelse(blood_type_abun_df$rs8176719 == "0", "T_T",
                                       ifelse(blood_type_abun_df$rs8176719 == "1", "TC_T",
                                              ifelse(blood_type_abun_df$rs8176719 == "2", "TC_TC", "err")))

blood_type_abun_df$rs8176746 <- ifelse(blood_type_abun_df$rs8176746 == "0", "G_G",
                                       ifelse(blood_type_abun_df$rs8176746 == "1", "T_G",
                                              ifelse(blood_type_abun_df$rs8176746 == "2", "T_T", "err")))

blood_type_abun_df$rs635634 <- ifelse(blood_type_abun_df$rs635634 == "0", "C_C",
                                      ifelse(blood_type_abun_df$rs635634 == "1", "T_C",
                                             ifelse(blood_type_abun_df$rs635634 == "2", "T_T", "err")))
blood_type_abun_df$rs7030248 <- ifelse(blood_type_abun_df$rs7030248 == "0", "G_G",
                                       ifelse(blood_type_abun_df$rs7030248 == "1", "A_G",
                                              ifelse(blood_type_abun_df$rs7030248 == "2", "A_A", "err")))
# blood_type_abun_df$rs1047781_t_t <- blood_type_abun_df$rs1047781
blood_type_abun_df$rs1047781_t <- ifelse(blood_type_abun_df$rs1047781_t == "0", "A_A",
                                         ifelse(blood_type_abun_df$rs1047781_t == "1", "T_A",
                                                ifelse(blood_type_abun_df$rs1047781_t == "2", "T_T", "err")))

# blood_type_abun_df$rs7030248
write.csv(blood_type_abun_df,file = 'blood_FUT2_2.csv',quote = F)

# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type.pdf", width = 12, height = 12)
# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type_FUT2all.pdf", width = 12, height = 12)
# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type_FUT2_non_secretors.pdf", width = 12, height = 12)
# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type_FUT2_secretors.pdf", width = 12, height = 12)
# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type_FUT2_secretors_AA.pdf", width = 12, height = 12)
# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type_FUT2_secretors_A_a.pdf", width = 12, height = 12)
# ggsave(filename = "blood_group_THSBC_all_newdata_p_genus_all_type_FUT2_non_secretors_aa.pdf", width = 12, height = 12)

# 计算 median ----
blood_type_abun_df %>% 
  group_by(rs1047781,Group) %>%
  summarise(across(starts_with("F.longum"):starts_with("g__Faecalibacterium.s__Fp"), median, na.rm = TRUE)) %>%
  ungroup() -> median_df1

blood_type_abun_df %>% 
  group_by(Group3) %>%
  summarise(across(starts_with("F.longum"):starts_with("g__Faecalibacterium.s__Fp"), median, na.rm = TRUE)) %>%
  ungroup() -> median_df2





if(norm!=T){
  save(blood_type_abun_df,file = '01_data/blood_type_abun_df.norm.F.RData')
}
if(norm!=T){
  save(blood_type_abun_df,file = '01_data/blood_type_abun_df.norm.F1.RData')
}
write.csv(blood_type_abun_df,file = 'blood_FUT2_1.csv',quote = F)

# library(ggplotify)
# P <- as.ggplot(P)
# 丰度 prevalence 展示 ----
feat_long <- abun[,-c(1)] %>%
  as.data.frame() %>%
  gather(key = "Species", value = "value")
feat_long$value <- as.numeric(feat_long$value)

feat_long_h <- feat_long %>% 
  # filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(value)
feat_long_h <- arrange(feat_long_h, desc(value))
feat_long_h_t <- feat_long
feat_long_h_t$Species <- factor(feat_long_h_t$Species, 
                                # levels = unique(feat_list_Faecalibacterium_df_long$Species[order(tapply(feat_list_Faecalibacterium_df_long$value, feat_list_Faecalibacterium_df_long$Species, median))])
                                levels = feat_long_h$Species
)

p3 <- ggboxplot(feat_long_h_t , x = 'Species', y ='value' ,
                # fill = 'Group',
                # fill = 'Group',
                # palette = pal_igv("default")(51)[-c(1,2)],
                width = 0.4,
                color =  "#386CB0",
                outlier.size = 0.5,
                
                outlier.colour = NA
)+
  
  # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
  # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
  scale_fill_manual(values =colors)+
  
  # scale_color_igv(alpha = 0.7)+
  # scale_fill_igv(alpha = 0.3)+
  # scale_y_log10()+
  # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
  #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
  stat_compare_means(
    ref.group = "F.longum",
    label = 'p.signif',
    method='wilcox.test',
    label.y = 3.8,
    paired = T
    # method.args = list(alternative = "less")
  )+
  # ylim(0,0.06)+
  # coord_cartesian(ylim = c(-0.001, 0.040)*100)+
  
  ylab("Abundance %")+xlab('')+
  
  # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
  theme_bw() +
  theme(
    legend.position="right",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(angle=90, hjust=1,face = 'italic',size=10),
    axis.text.y=element_text(face = 'plain',size=10),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )


# p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p3)


ggsave(filename = "F.longum.box.BGI.pdf", width = 4.5, height = 5)
# ggsave(filename = "F.box.non.pdf", width = 5, height = 7)

library(tidyr)
#abundance
feat_long_a <- feat_long %>% 
  # filter(Group == 'Control') %>%
  group_by(Species) %>%
  summarise(abundance = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(abundance))
#prevalence
feat_long_p <- feat_long %>% 
  group_by(Species) %>%
  summarise(prevalence = sum(value > 0, na.rm = TRUE) / n()) %>%
  ungroup() %>%
  arrange(desc(prevalence))

feat_long_ap <- merge(feat_long_a, feat_long_p, by = "Species")
feat_long_ap <- arrange(feat_long_ap, desc(abundance))
feat_long_ap$Species <- factor(feat_long_ap$Species, levels = feat_long_ap$Species)
feat_long_ap$Species

pA <- ggplot(feat_long_ap, aes(x = Species,  y = abundance)) + #, fill = bacterium
  # geom_bar(stat = "identity",fill='#B3E2CD') + #reorder(bacterium,-mean_abundance)
  geom_bar(stat = "identity",fill='#B3E2CD') + #reorder(bacterium,-mean_abundance)
  ylab("Abundance (mean) %") +
  xlab("") +
  # geom_hline(yintercept = 0.5,color='red')+
  theme(legend.position = 'none')+
  theme_bw() +
  theme(
    legend.position="none",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.y=element_text(face = 'plain',size=10),
    axis.text.x=element_text(angle=90,vjust = 0.5, hjust=1,face = 'italic',size=10),
    # axis.text.x=element_blank(),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )
pA
pP <- ggplot(feat_long_ap, aes(x = Species, y = prevalence*100)) + #fill = bacterium
  # geom_bar(stat = "identity", fill = "#FDCDAC") + #x = reorder(bacterium,-prevalence)
  geom_bar(stat = "identity", fill = "#C9E3FF") + #x = reorder(bacterium,-prevalence)
  ylab("Prevalence %") +
  xlab("") +
  # geom_hline(yintercept = 25,color='red')+
  # ggtitle("Prevalence of Bacteria") +
  theme(legend.position = 'none')+
  theme_bw() +
  theme(
    legend.position="none",
    # axis.text.x=element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.y=element_text(face = 'plain',size=10),
    axis.text.x=element_blank(),
    text = element_text(size=12,face = 'plain',family ='',colour = 'black')
  )
pP

p <- ggarrange(pP,pA, ncol = 1, nrow = 2,
               # labels = c("Prevalence", "B"),
               common.legend = TRUE, legend = "right",
               heights  = c(1, 1.6))
p
ggsave(filename = "F.longum.prevalence.abundance.BGI.pdf", plot = p, width = 5, height = 4.5)





# 合并data ----
a <- intersect(blood_type_df1$sample,blood_type_df2$sample)
blood_type_abun_df1_t <- blood_type_abun_df1[! (blood_type_abun_df1$sample %in% a),]
blood_type_abun_df1_t <- blood_type_abun_df1_t[,-c(2:3)]  # 去掉样本信息
# blood_type_abun_df2
blood_type_abun_df2_t <- blood_type_abun_df2[,-c(2:4)]
colnames(blood_type_abun_df1_t)
colnames(blood_type_abun_df2_t)
blood_type_abun_df12_t <- rbind(blood_type_abun_df1_t, blood_type_abun_df2_t)

p <- list()
for (species in unique(colnames(blood_type_abun_df12_t)[-c(1:4)])) {
  p[[species]] <- ggboxplot(blood_type_abun_df12_t , x = 'Group', y =species ,
                            fill = 'rs1047781',
                            # palette = pal_igv("default")(51)[-c(1,2)],
                            width = 0.4,
                            # color =  "#386CB0",
                            outlier.size = 0.5,
                            
                            outlier.colour = NA
  )+
    facet_grid(~rs1047781, scales = 'free') +
    # rs1047781
    # scale_color_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,2,4,9)])+
    # scale_fill_manual(values =pal_igv("default",alpha = 0.8)(51)[-c(1,2,4,9)])+
    # scale_fill_manual(values =pal_igv("default",alpha = 1)(51)[-c(1,4,9)])+
    scale_fill_manual(values =colors)+
    # stat_pvalue_manual(stat.test,label = "p")+
    # scale_color_igv(alpha = 0.7)+
    # scale_fill_igv(alpha = 0.3)+
    # scale_y_log10()+
    # stat_summary(fun = "median", geom = "text", aes(label = round(..y.., 2)),
    #              position = position_dodge(width = 0.25), vjust = 2) +  # 不再需要指定 y 值
    # stat_compare_means(
    #   ref.group = "F.longum",
    #                    label = 'p.signif',
    #                    method='wilcox.test',
    #                    label.y = 4.7,
    #                    paired = T
    #                    # method.args = list(alternative = "less")
    # )+
    stat_compare_means( 
      label = 'p',# 'p.signif'
      #                    method='wilcox.test',
      # label.y = 6,
      method='wilcox.test'
      # method.args = list(alternative = "less") #  greater
    )+
    # ylim(0,0.06)+
    # coord_cartesian(ylim = c(-0.001, 0.08)*100)+
    
    ylab("Abundance %")+xlab('')+
    
    # scale_y_continuous(labels = trans_format("log10", math_format(10^.x))) +  # 对 y 轴进行10的对数变换
    theme_bw() +
    theme(
      legend.position="none",
      # axis.text.x=element_blank(),
      # axis.ticks.x = element_blank(),
      axis.text.x=element_text( hjust=1,face = 'plain',size=10),
      axis.text.y=element_text(face = 'plain',size=10),
      text = element_text(size=12,face = 'plain',family ='',colour = 'black')
    )+
    ggtitle(species)
}
P <- ggarrange(p$F.longum, 
               p$F.prausnitzii,
               p$F.faecis_F.taiwanense, 
               p$F.duncaniae, 
               p$F.hattorii,
               p$t__SGB15323,
               p$t__SGB15317,
               p$t__SGB15326,
               p$t__SGB15333,
               # p$g__Faecalibacterium,
               # p$g__Faecalibacterium.s__Fp,
               ncol = 3,
               nrow = 3)

P <- ggarrange(p$F.longum, 
               p$F.prausnitzii,
               p$F.faecis_F.taiwanense, 
               p$F.duncaniae, 
               p$F.hattorii,
               p$t__SGB15323,
               p$t__SGB15317,
               p$t__SGB15326,
               p$t__SGB15333,
               p$g__Faecalibacterium,
               p$g__Faecalibacterium.s__Fp,
               ncol = 4,
               nrow = 3)
print(P)


# 计算median ----
blood_type_abun_df1 %>% 
  group_by(rs1047781,Group2) %>%
  summarise(across(starts_with("F.prausnitzii"):starts_with("g__Faecalibacterium.s__Fp"), median, na.rm = TRUE)) %>%
  ungroup() -> median_df1

blood_type_abun_df1 %>% 
  group_by(rs1047781,Group1) %>%
  summarise(across(starts_with("F.prausnitzii"):starts_with("g__Faecalibacterium.s__Fp"), mean, na.rm = TRUE)) %>%
  ungroup() -> mean_df1

blood_type_abun_df1 %>% 
  group_by(Group2) %>%
  summarise(across(starts_with("F.prausnitzii"):starts_with("g__Faecalibacterium.s__Fp"), median, na.rm = TRUE)) %>%
  ungroup() -> median_df1

blood_type_abun_df1 %>% 
  group_by(Group2) %>%
  summarise(across(starts_with("F.prausnitzii"):starts_with("g__Faecalibacterium.s__Fp"), mean, na.rm = TRUE)) %>%
  ungroup() -> mean_df1

blood_type_abun_df2_t %>% 
  group_by(rs1047781) %>%
  summarise(across(starts_with("F.prausnitzii"):starts_with("g__Faecalibacterium.s__Fp"), median, na.rm = TRUE)) %>%
  ungroup() -> median_df2







# Recall Precision ----
library(readxl)
PR <- read_excel("./01_data/P_R.xlsx",sheet = 1)
colnames(PR) <- c('Species','Precision','Recall','all')
y_label <- round(PR$all[1],2)*100
PR$all <- NULL
PR$Precision <- round(PR$Precision,2)*100
PR$Recall <- round(PR$Recall,2)*100

PR <- PR %>% gather(Category,Value,-Species)

# bar 图 
PR$Species <- factor(PR$Species,levels=c("F.prausnitzii","F.duncaniae","F.longum","F.others"))
ggplot(PR, aes(x = Species, y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7,position = "dodge") +
  theme_minimal() +
  labs(title = "", x = "", y = "Value", fill = "Group") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # 标注每组的总和
  geom_text(
    aes(label = Value, group = Species,y= Value/2),
    position = position_dodge2(width = 0.9),
    size = 5,
    color = "white"
  ) +
  # 标注所有分组中的最大值
  # geom_hline(yintercept = max(tapply(sum_df$sum_62, sum_df$Group, sum)), color = "red") +
  # annotate(
  #   "text",
  #   x = 1,
  #   y = max(tapply(sum_df$sum_62, sum_df$Group, sum)) + 0.5,
  #   label = paste("Max: ", max(tapply(sum_df$sum_62, sum_df$Group, sum))),
  #   color = "red"
  # )+
  # 设置自定义颜色方案
  scale_fill_manual(
    values = c(
      # "#FDBF6F", "#FF7F00"
      # "#A6CEE3", "#1F78B4"
      "#33A02C","#B2DF8A" 
      # "#A6CEE3", "#1F78B4"
    )
  )+
  ylim(0,120)+
  # scale_fill_manual(
  #   values = c(
  #     "#56B4E9", "#009E73", "#F0E442",
  #     "#D55E00", "#CC79A7", "#E69F00",
  #     "#0072B2", "#999999"
  #   )
  # )+
  # scale_y_continuous(breaks = seq(4.5,6,0.5))+
  geom_hline(yintercept = y_label, color = "red", linetype = "dashed", size = 0.5) +
  # 为横线添加文字标注
  annotate("text", 
           x = Inf,  # 靠右显示（Inf表示最右侧，-Inf表示最左侧）
           y = y_label+10,  # 与横线y值一致
           label = paste0('Overall Accuracy: ', y_label,'%'),  # 标注文字
           color = "black", 
           hjust = 1.1,  # 水平调整（1.1表示在横线右侧）
           vjust = -0.5,  # 垂直调整（-0.5表示在横线上方）
           size = 4) +
  
  theme_bw()+
  theme(legend.position = 'right')+
  theme(
    text = element_text(family = "",size = 12,face = 'plain'),
    # axis.line=element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text( hjust=0.5, vjust=0,family = '',size = 10,face = 'italic'),#angle=90,
    # axis.text.x =element_blank(),
    # axis.ticks.x = element_blank(),
    # panel.grid=element_blank(),
    # panel.border=element_blank(),
    # strip.background = element_blank(),
    # strip.text = element_blank(),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 12),
    # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),     legend.position = 'right'
  )
ggsave('PR.pdf',width=8,height = 6)


















