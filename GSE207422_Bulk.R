rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE207422_NSCLC/")

library(data.table)
library(tidyverse)
library(dplyr)
library(ggpubr)

pdata <- read.table("GSE207422_NSCLC_bulk_RNAseq_metadata.csv", header = T, sep = ",")

pdata$RECIST <- factor(pdata$RECIST, levels = c("SD", "PR", "CR"))
pdata <- pdata %>% arrange(RECIST)

tpm <- read.table("GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", header = T, sep = "\t")
tpm <- column_to_rownames(.data = tpm, var = "Gene")

tpm <- tpm[,pdata$Sample]
pdata$Sample
colnames(tpm)

# SD:7
# PR:14
# CR:3

group_list = c(rep("SD", 7), rep("PR", 14), rep("CR", 3))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "SD")

library(limma) 
data = normalizeBetweenArrays(tpm)
#去除批次效应后
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-1,4))
#去除批次效应前data#去除批次效应前
boxplot(tpm, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,6))

###箱线图
data <- as.data.frame(data)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("SD", 7), rep("PR", 14), rep("CR", 3))

z$type <- factor(z$type, levels = c("SD", "PR", "CR"))

colnames(z)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") +
  ylab(paste("Gene Expression log2 TPM")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("SD", "PR"), c("SD", "CR")), 
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(4.5, 5)
  ) + ggtitle("IGSF3")


ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)
