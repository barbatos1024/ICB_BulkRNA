rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE302922_CRC/")

library(data.table)
library(tidyverse)
library(dplyr)
library(ggpubr)

tissue <- read.table("GSE302922_Tissue_count_coding.csv.gz", header = T, sep = ",")
tissue <- distinct(tissue, Gene.name, .keep_all = T)
rownames(tissue) <- tissue$Gene.name
tissue <- tissue[,-c(1,2,3)]
colnames(tissue)

# Pre_R: 1,3,5,7,10,12
# Pre_NR:  2,4,6,8,9,11

tissue <- tissue[,c(1,3,5,7,10,12, 2,4,6,8,9,11)]   # 排序 Pre_R, Pre_NR

### 创建分组 ###
group <- c(rep("Pre_R", 6), rep("Pre_NR", 6)) 

library(edgeR)

### 数据预处理
# 构建 DGEList 对象
dgelist <- DGEList(counts = tissue, group = group)
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
z$type = c(rep("Pre_R", 6), rep("Pre_NR", 6)) 

z$type <- factor(z$type, levels = c("Pre_R", "Pre_NR"))

colnames(z)

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
    label.y = c(6.8)
  ) + ggtitle("IGSF3")

ggsave("IGSF3_tissue.png", width = 3, height = 5, dpi = 300)

########################################################
####################################################################################
########################################################
####################################################################################
########################################################

evrna <- read.table("GSE302921_evRNA_count_coding.csv.gz", header = T, sep = ",")
evrna <- distinct(evrna, Gene.name, .keep_all = T)
rownames(evrna) <- evrna$Gene.name
evrna <- evrna[,-c(1,2,3)]
colnames(evrna)

pdata <- read.table("evRNA_characteristic.csv", header = T, sep = ",")
pdata$Condition <- paste(pdata$Treatment, pdata$Response, sep = "_")
pdata <- arrange(.data = pdata, Condition)

colnames(evrna) %in% pdata$Sample

evrna <- evrna[,pdata$Sample]     # 根据pdata的样本顺序进行排序

table(pdata$Condition)
### 创建分组 ###
group <- c(rep("Post_PR", 12), rep("Post_nonPR", 13), rep("Pre_PR", 12), rep("Pre_nonPR", 13)) 

library(edgeR)

### 数据预处理
# 构建 DGEList 对象
dgelist <- DGEList(counts = evrna, group = group)
# 过滤 low count 数据 例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1) >= 0                ###要求至少在2个样本中表达量大于1
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
z$type = c(rep("Post_PR", 12), rep("Post_nonPR", 13), rep("Pre_PR", 12), rep("Pre_nonPR", 13)) 

z$type <- factor(z$type, levels = c("Pre_PR", "Post_PR", "Pre_nonPR", "Post_nonPR"))

colnames(z)

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
    comparisons = list(c("Pre_PR", "Post_PR"), 
                       c("Pre_nonPR", "Post_nonPR"),
                       c("Post_PR", "Post_nonPR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(1,1,1.5)
  ) + ggtitle("IGSF3")

ggsave("IGSF3_evrna.png", width = 3, height = 5, dpi = 300)



