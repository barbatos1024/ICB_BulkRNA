rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE126044_NSCLC/")

library(data.table)
library(tidyverse)
library(dplyr)

counts <- read.table("GSE126044_counts.txt", header = T, sep = "\t")
counts <- distinct(counts, X, .keep_all = T)
rownames(counts) <- counts$X
counts <- counts[,-1]
colnames(counts)

# Pre_R: "Dis_02", "Dis_15", "Dis_04", "Dis_17", "Dis_10"
# Pre_NR: 



counts <- counts[,c(2,5,6,8,11, 1,3,4,7,9,10,12,13,14,15,16)]   # 排序 Pre_R, Pre_NR

### 创建分组 ###
group <- c(rep("Pre_R", 5), rep("Pre_NR", 11)) 

library(edgeR)

### 数据预处理
# 构建 DGEList 对象
dgelist <- DGEList(counts = counts, group = group)
# 过滤 low count 数据 例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1) >= 2                ###要求至少在2个样本中表达量大于1
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]                  
# 标准化 以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

# 获取 log2 CPM 表达矩阵
log_cpm <- cpm(dgelist_norm, log = TRUE)

# 画箱型图
boxplot(log_cpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 CPM Boxplot", ylab = "Log2 CPM")        ########！！！！！！！


###箱线图
logdata <- as.data.frame(log_cpm)
GENE = logdata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 5), rep("Pre_NR", 11)) 

z$type <- factor(z$type, levels = c("Pre_R", "Pre_NR"))

colnames(z)

library(ggpubr)


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log cpm")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Pre_R", "Pre_NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7)
  ) + ggtitle("IGSF3")

ggsave("ISGF3_PD-1.png", width = 3, height = 5, dpi = 300)

