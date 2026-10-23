
rm(list = ls())

library(data.table)
library(ggplot2)
library(ggpubr)
library(org.Hs.eg.db)
library(tidyverse)
library(ggpubr)

setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/PRJEB25780_PD1_Gastric/")

data <- read.table("ICB.Kim2018_Pembrolizumab_Gastric.self_subtract.txt", 
                   header = T, check.names = F)

entrez_ids <- as.character(data$Entrez)

# 转换
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

gene_symbols
data$genesymbol <- gene_symbols

data <- distinct(.data = data, genesymbol, .keep_all = T)
data <- na.omit(data)
rownames(data) <- data$genesymbol
colnames(data)
data <- data[,c(2:46)]


pdata <- read.table("ICB.Kim2018_Pembrolizumab_Gastric.clinical.txt", 
                    header = T, sep = "\t")

pdata_filtered <- pdata %>%
  filter(patient %in% colnames(data))

colnames(data) %in% pdata_filtered$patient

a = pdata_filtered$patient

data <- data[,a]

group_list = c(rep('R', 12), rep("NR", 33))
group_list = factor(group_list)
group_list <- relevel(group_list, ref = "R")

library(limma)

##除去批次效应
nordata = normalizeBetweenArrays(data)   # 已经 log2(data+1) 处理
#去除批次效应后
boxplot(nordata, outline = FALSE, notch = T, col = group_list, las = 2)
#去除批次效应前
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2)

### 箱线图
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('R', 12), rep("NR", 33))

z$type <- factor(z$type, levels = c("R", "NR"))

colnames(z)


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Relative Gene Expression")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c('R', 'NR')), 
    label.y = c(2.2), method = "t.test"
  ) + ggtitle("IGSF3")

ggsave("IGSF3.PNG", width = 3, height = 5, dpi = 300)


