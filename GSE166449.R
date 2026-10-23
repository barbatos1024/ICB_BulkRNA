rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE166449_NSCLC/")

library(data.table)
library(tidyverse)
library(dplyr)

tpm <- read.table("GSE166449_Raw_gene_TPM_matrix.txt", header = T, sep = "\t", row.names = 1)
colnames(tpm)

### 创建分组 ###
group <- c(rep("Pre_R", 7), rep("Pre_NR", 15)) 

library(limma)
logtpm <- log(tpm + 1)

nortpm <- normalizeBetweenArrays(logtpm)

boxplot(nortpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)") 

boxplot(logtpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)")     


###箱线图
nortpm <- as.data.frame(nortpm)
GENE = nortpm["IGSF3",]
z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 7), rep("Pre_NR", 15)) 

z$type <- factor(z$type, levels = c("Pre_R", "Pre_NR"))

colnames(z)

library(ggpubr)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2(TPM + 1)")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Pre_R", "Pre_NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(1.5)
  ) + ggtitle("IGSF3")

ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)
