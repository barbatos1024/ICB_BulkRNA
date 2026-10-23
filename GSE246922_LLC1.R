
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE246922_LLC1_KP/")

library(data.table)
library(tidyverse)
library(dplyr)


llc1 <- read.table("GSE246922_LLC1_RNA_counts_vst.csv", header = T, sep = ",")
colnames(llc1)[1] <- "ensembl_gene_id"

library(biomaRt)
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_ids <- llc1$ensembl_gene_id

annotations <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mouse
)

merged_df <- merge(llc1, annotations, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(merged_df)

merged_df <- distinct(merged_df, mgi_symbol, .keep_all = T)
merged_df <- na.omit(merged_df)

rownames(merged_df) <- merged_df$mgi_symbol
colnames(merged_df)

data <- merged_df[,c(2:10)]
colnames(data)


### 创建分组 ###
group <- c(rep("LLC1_Parental", 3), rep("LLC1_Ifng", 3), rep("LLC1_Relapse", 3))

library(edgeR)

### 数据预处理
# 构建 DGEList 对象
dgelist <- DGEList(counts = data, group = group)
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
nordata <- as.data.frame(log_cpm)
GENE = nordata["Igsf3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("LLC1_Parental", 3), rep("LLC1_Ifng", 3), rep("LLC1_Relapse", 3))

z$type <- factor(z$type, levels = c("LLC1_Parental", "LLC1_Ifng", "LLC1_Relapse"))

colnames(z)

library(ggpubr)

ggplot(z, aes(type, Igsf3, fill = type)) +
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
    comparisons = list(c("LLC1_Parental", "LLC1_Ifng"),
                       c("LLC1_Parental", "LLC1_Relapse")),
    label.y = c(5.32, 5.35), method = "t.test"  ### 或用"t.test" 进行配对非参数检验
  ) + ggtitle("Igsf3")


