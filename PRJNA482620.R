
rm(list = ls())

getwd()
setwd('/Users/naturexv/Datasets/PD-1/BulkRNA/PRJNA482620_GBM/')

library(data.table)
library(limma)
library(tidyverse)
library(ggpubr)

expr <- read.table("GBM-PRJNA482620.expression_matrix.tsv", header = T)
rownames(expr) <- expr$GENE_SYMBOL
expr <- expr[,-1]

colnames(expr) %in% pdata$sample_id

pdata <- read.table("GBM-PRJNA482620.clinical_characteristic.tsv", header = T, sep = "\t")
pdata <- arrange(.data = pdata, response)

expr <- expr[,pdata$sample_id]   # 将expr样本排列顺序按照pdata排列（N,R）
table(pdata$response)

group_list <- c(rep('NR', 17), rep('R', 17))
group_list = factor(group_list)
group_list <- relevel(group_list, ref = "R")

data <- log2(expr + 1)

##除去批次效应
data_normalized = normalizeBetweenArrays(data)   # 已经 log2(data+1) 处理
#去除批次效应后
boxplot(data_normalized, outline = FALSE, notch = T, col = group_list, las = 2)
#去除批次效应前
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2)

### 箱线图
data_normalized <- as.data.frame(data_normalized)
GENE = data_normalized["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('pre_NR', 17), rep('pre_R', 17))

z$type <- factor(z$type, levels = c("pre_R", "pre_NR"))

colnames(z)


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2(RPKM + 1)")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c('pre_R', 'pre_NR')), 
    label.y = c(6), method = "t.test"
  ) + ggtitle("IGSF3")

ggsave("IGSF3.PNG", width = 4, height = 5, dpi = 300)






