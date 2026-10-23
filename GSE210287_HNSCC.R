rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE210287_HNSCC/")

library(data.table)
library(tidyverse)
library(dplyr)

counts <- read.table("GSE210287_raw_counts_18_0606.csv", header = T, sep = ",", row.names = 1)
colnames(counts)
# Pre_NR: "S5"  "S7"  "S11" "S13" "S22" 
# Pre_R: "S9"  "S16" "S20" "S70" "S21" "S24" "S29"
# Post_NR: "S6"  "S8"  "S12" "S23" 
# Post_R: "S10" "S15" "S17" "S19" "S71" "S25"

counts <- counts[,c(6:12, 1:5, 17:22, 13:16)]   # 排序 Pre_R, Pre_NR, Post_R, Post_NR

### 创建分组 ###
group <- c(rep("Pre_R", 7), rep("Pre_NR", 5), rep("Post_R", 6), rep("Post_NR", 4)) 

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
z$type = c(rep("Pre_R", 7), rep("Pre_NR", 5), rep("Post_R", 6), rep("Post_NR", 4)) 

z$type <- factor(z$type, levels = c("Pre_R","Post_R", "Pre_NR", "Post_NR"))

colnames(z)

library(ggpubr)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") +
  ylab(paste("Gene Expression log2 (FPKM + 1)")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("Pre_R","Post_R"), c("Pre_NR", "Post_NR")), 
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7, 7)
  ) + ggtitle("IGSF3")


ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)
