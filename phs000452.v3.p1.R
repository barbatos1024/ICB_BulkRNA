
rm(list = ls())

library(data.table)
library(ggplot2)
library(ggpubr)
library(org.Hs.eg.db)
library(tidyverse)
library(ggpubr)

setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/phs000452.v3.p1/")

data <- read.table("phs000452.v3.p1_RNA-seq TPM matrix.txt", header = T, check.names = F)
data <- as.data.frame(t(data))

pdata <- read.table("phs000452.v3.p1_clinical_information.tsv", 
                    header = T, sep = "\t")

pdata_filtered <- pdata %>%
  filter(Sample.ID %in% colnames(data))

colnames(data) %in% pdata_filtered$Sample.ID
table(pdata_filtered$BR)
# CR 16
# MR 2
# PD 56
# PR 31
# SD 16
table(pdata_filtered$priorCTLA4)
# naive 74
# post 47


# 根据response排序
pdata_filtered <- arrange(pdata_filtered, BR)
a = pdata_filtered$Sample.ID

data <- data[,a]     # 根据response排序

group_list = c(rep('CR', 16), rep("MR/SD", 2), rep("PD", 56), rep("PR", 31), rep("MR/SD", 16))
group_list = factor(group_list)
group_list <- relevel(group_list, ref = "CR")

library(limma)

logdata <- log2(data + 1)

##除去批次效应
nordata = normalizeBetweenArrays(logdata)  
#去除批次效应后
boxplot(nordata, outline = FALSE, notch = T, col = group_list, las = 2)
#去除批次效应前
boxplot(logdata, outline = FALSE, notch = T, col = group_list, las = 2)

### 箱线图
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
# z$type = c(rep('CR', 16), rep("MR/SD", 2), rep("PD", 56), rep("PR", 31), rep("MR/SD", 16))
# z$type <- factor(z$type, levels = c("CR", "PR", "MR/SD", "PD"))

z$type = c(rep('Responder', 16), rep("Responder", 2), rep("Progressor", 56), rep("Responder", 31), rep("Responder", 16))
z$type <- factor(z$type, levels = c("Responder", "Progressor"))

colnames(z)


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 (TPM + 1)")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Responder", "Progressor")), 
    label.y = c(8.4), method = "t.test"
  ) + ggtitle("IGSF3")

ggsave("IGSF3_pre_Respond_Progress.PNG", width = 3, height = 5, dpi = 300)



############# 根据ctla4治疗排序 #############
pdata_filtered <- arrange(pdata_filtered, priorCTLA4)
a = pdata_filtered$Sample.ID

data <- data[,a]     # 根据response排序

group_list = c(rep('CTLA4naive', 74), rep("CTLA4treated", 47))
group_list = factor(group_list)
group_list <- relevel(group_list, ref = "CTLA4naive")

library(limma)

logdata <- log2(data + 1)

##除去批次效应
nordata = normalizeBetweenArrays(logdata)  
#去除批次效应后
boxplot(nordata, outline = FALSE, notch = T, col = group_list, las = 2)
#去除批次效应前
boxplot(logdata, outline = FALSE, notch = T, col = group_list, las = 2)

### 箱线图
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框

z$type = c(rep('CTLA4naive', 74), rep("CTLA4treated", 47))
z$type <- factor(z$type, levels = c("CTLA4naive", "CTLA4treated"))

colnames(z)


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 (TPM + 1)")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("CTLA4naive", "CTLA4treated")), 
    label.y = c(8.4), method = "t.test"
  ) + ggtitle("IGSF3")

ggsave("IGSF3_CTLA4naive_CTLA4treated.PNG", width = 3, height = 5, dpi = 300)










