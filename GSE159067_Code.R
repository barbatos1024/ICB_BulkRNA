
rm(list = ls())   # 无IGSF3基因

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE159067_HNSCC/")

library(data.table)
library(tidyverse)

library(edgeR)

pdata <- fread('pdata.csv', header = T)
pdata$patient
pdata$response
pdata$rank <- c(paste("s", c(1:102), sep = ""))

log_cpm <- read.table('GSE159067_IHN_log2cpm_data.txt.gz', header = T, row.names = 1)
colnames(log_cpm) <- pdata$rank

pdata <- arrange(.data = pdata, response)
rank <- pdata$rank

log_cpm <- log_cpm[,c(rank)]   # 按照response进行排序
table(pdata$response)


### 创建分组 ###
group <- c(rep("CR",5), rep("PD", 64), rep("PR", 6), rep("SD", 27))


# 画箱型图
boxplot(log_cpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 CPM Boxplot", ylab = "Log2 CPM")        ########！！！！！！！


###箱线图
nordata <- as.data.frame(log_cpm)
GENE = nordata["CD274",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("CR",5), rep("PD", 64), rep("PR", 6), rep("SD", 27))

z$type <- factor(z$type, levels = c("CR", "PR", "SD", "PD"))

colnames(z)

library(ggpubr)

ggplot(z, aes(type, CD274, fill = type)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log cpm")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("HC", "LN")),
    label.y = c(7.5), method = "t.test"  ### 或用"t.test" 进行配对非参数检验
  ) + ggtitle("CD274")


