rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE96619_SKCM/")

library(data.table)
library(tidyverse)
library(ggrepel)
library(ggpubr)

fpkm <- read.table("GSE96619_PatientFPKM.csv", header = T, sep = ",")
fpkm <- distinct(fpkm, Gene, .keep_all = T)
fpkm <- column_to_rownames(fpkm, var = "Gene")
colnames(fpkm)

data <- fpkm[,c(2,4,1,3,6,8,10,5,7,9)]

# 创建分组
group_list <- factor(c(rep("Pre_R", 2), rep("Post_R", 2), rep("Pre_NR", 3), rep("Post_NR", 3)))
group_list <- relevel(group_list, ref = "Pre_R")

library(limma)
logdata <- log2(data + 1)
nordata <- normalizeBetweenArrays(logdata)
#去除批次效应后
boxplot(nordata, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,10))

#去除批次效应前
boxplot(logdata, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,10))

###箱线图
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 2), rep("Post_R", 2), rep("Pre_NR", 3), rep("Post_NR", 3))

z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))

colnames(z)

### 加个p值 ####### 配对 #######
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 2), rep("Post_R", 2), rep("Pre_NR", 3), rep("Post_NR", 3))

z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))

z$pairinfo = c(rep(1:2, 2), rep(3:5, 2))


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5, colour = "red") +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 (FPKM + 1)")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("Pre_R","Post_R"), c("Pre_NR", "Post_NR")),
    paired = TRUE, method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7, 6.5)
  ) + ggtitle("IGSF3")
###   geom_jitter(shape = 16, position = position_jitter(0.2)) + #散点抖动 ###

ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)




