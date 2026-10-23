
rm(list = ls())

setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/MGH_SKCM/")
getwd()

list.files()

library(data.table)
library(tidyverse)


### GSE115821
GSE115821 <- read.table("GSE115821_MGH_counts.csv", header = T, sep = ",")
GSE115821 <- distinct(.data = GSE115821, Geneid, .keep_all = T)
rownames(GSE115821) <- GSE115821$Geneid
GSE115821 <- GSE115821[,-c(1:6)]
colnames(GSE115821)

GSE115821_pdata <- read.table("GSE115821_pdata.csv", header = T, sep = ",")
GSE115821_pdata$Sample

GSE115821_pdata$Sample %in% colnames(GSE115821)   # 确认两样本名是否完全重合
GSE115821_pdata$Batch <- c(rep("GSE115821", 37))
GSE115821_pdata$GSE <- c(rep("GSE115821", 37))
### GSE115821
GSE168204 <- read.table("GSE168204_MGH_counts.csv", header = T, sep = ",", check.names = F)
GSE168204 <- distinct(.data = GSE168204, Gene, .keep_all = T)
GSE168204 <- column_to_rownames(.data = GSE168204, var = "Gene")

GSE168204_pdata <- read.table("GSE168204_pdata.csv", header = T, sep = ",")

GSE168204_pdata$Sample %in% colnames(GSE168204)   # 确认两样本名是否完全重合
GSE168204_pdata$Batch <- c(rep("GSE168204", 27))
GSE168204_pdata$GSE <- c(rep("GSE168204", 27))
### GSE273583
GSE273583 <- read.table("GSE273583_counts.csv", header = T, sep = ",", check.names = F)
GSE273583 <- distinct(.data = GSE273583, Geneid, .keep_all = T)
GSE273583 <- column_to_rownames(.data = GSE273583, var = "Geneid")

GSE273583_pdata <- read.table("GSE273583_pdata.csv", header = T, sep = ",")
GSE273583_pdata$GSE <- c(rep("GSE273583", 16))

### pdata 整理
pdata_list <- list(GSE273583_pdata, GSE168204_pdata, GSE115821_pdata)
metadata <- rbindlist(pdata_list, use.names = TRUE, fill = TRUE)

metadata$Treatment <- gsub("PRE", "Pre", metadata$Treatment, ignore.case = TRUE)
metadata$Treatment <- gsub("POST", "Post", metadata$Treatment, ignore.case = TRUE)
metadata$Treatment <- trimws(metadata$Treatment) # 去掉两端空格

### count, metadata 合并
metadata$Sample_GSE <- paste0(metadata$GSE, "_", metadata$Sample)

colnames(GSE115821) <- paste0("GSE115821_", colnames(GSE115821))
colnames(GSE168204) <- paste0("GSE168204_", colnames(GSE168204))
colnames(GSE273583) <- paste0("GSE273583_", colnames(GSE273583))

common_genes <- intersect(rownames(GSE115821), 
                          intersect(rownames(GSE168204), rownames(GSE273583)))

GSE115821 <- GSE115821[common_genes, ]
GSE168204 <- GSE168204[common_genes, ]
GSE273583 <- GSE273583[common_genes, ]


combined_counts <- cbind(GSE115821, GSE168204, GSE273583)
###

metadata$condition <- paste(metadata$Treatment, metadata$Response, sep = "_")
table(metadata$condition)
metadata$condition <- sub("^On", "Post", metadata$condition, ignore.case = TRUE)
table(metadata$condition)

table(metadata$Batch)

expr_log2 <- log2(combined_counts + 1)

# 2. 强制转为数值矩阵并清理非法值 (针对之前报错的防御)
final_mat <- as.matrix(expr_log2)
mode(final_mat) <- "numeric"
final_mat[!is.finite(final_mat)] <- 0

row_vars <- apply(final_mat, 1, var)
final_mat <- final_mat[row_vars > 0, ]

pca <- prcomp(t(final_mat), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Batch = metadata$Batch,
  Condition = metadata$condition
)

ggplot(pca_df, aes(PC1, PC2, color = Batch, shape = Condition)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCA Analysis: Checking Batch Effects",
       x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"))



library(sva)     # 去除批次效应
# count_matrix: 你的 combined_counts (必须是纯数字矩阵)
# batch: metadata$Batch
# group: metadata$condition
counts_mat <- as.matrix(combined_counts)
counts_mat[is.na(counts_mat)] <- 0
adjusted_counts <- ComBat_seq(counts = as.matrix(counts_mat), 
                              batch = metadata$GSE, 
                              group = metadata$condition)
# 对矫正后的数据做 Log2 转化
adj_expr_log2 <- log2(adjusted_counts + 1)

# 过滤方差并运行 PCA (复用之前的逻辑)
adj_final_mat <- as.matrix(adj_expr_log2[apply(adj_expr_log2, 1, var) > 0, ])
adj_pca <- prcomp(t(adj_final_mat), scale. = TRUE)

# 绘图查看
adj_pca_df <- data.frame(PC1 = adj_pca$x[,1], PC2 = adj_pca$x[,2], 
                         Batch = metadata$GSE)
ggplot(adj_pca_df, aes(PC1, PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA After ComBat-seq Correction")


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = round(adjusted_counts), # 确保是整数
                              colData = metadata,
                              design = ~ condition)

dds <- DESeq(dds)

# 提取归一化后的数据
vsd <- vst(dds, blind = FALSE)
plot_data <- plotCounts(dds, gene = "IGSF", intgroup = "condition", returnData = TRUE)

ggplot(plot_data, aes(x = condition, y = count, fill = condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_bw() +
  scale_y_log10() +
  labs(title = "Corrected IGSF3 Expression", y = "Normalized Counts")










