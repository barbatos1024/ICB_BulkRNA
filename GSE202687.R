rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE202687_SKCM/")

library(data.table)
library(tidyverse)
library(dplyr)

counts <- read.table("GSE202687_raw_counts_GRCh38.p13_NCBI.tsv", header = T, sep = "\t")

library(org.Hs.eg.db)
library(dplyr)

geneid <- counts$GeneID

mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = as.character(geneid),
                                 columns = c("SYMBOL"),
                                 keytype = "ENTREZID")
counts$genesymbol <- mapping$SYMBOL
colnames(counts)

counts <- distinct(.data = counts, genesymbol, .keep_all = T)
counts <- na.omit(counts)
rownames(counts) <- counts$genesymbol
counts <- counts[,-c(1,17)]

counts <- counts[,-c(3, 8, 13)]   # 去除未配对样本
counts <- counts[,c(1,3,5,7,9,11, 2,4,6,8,10,12)]   # pre, post

### 创建分组 ###
group <- c(rep("Pre", 6), rep("Post", 6)) 

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
GENE = logdata["IGSF3",]   # CD274 差异明显

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre", 6), rep("Post", 6)) 

z$type <- factor(z$type, levels = c("Pre", "Post"))
z$pairinfo <- rep(1:6, 2)

colnames(z)

library(ggpubr)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5, colour = "red") +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 cpm")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("Pre", "Post")), 
    paired = TRUE, method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(8)
  ) + ggtitle("IGSF3")
###   geom_jitter(shape = 16, position = position_jitter(0.2)) + #散点抖动 ###

ggsave("IGSF3.png", width = 4, height = 5, dpi = 300)
