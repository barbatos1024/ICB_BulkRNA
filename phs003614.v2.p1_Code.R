
rm(list = ls())

setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/phs003614.v2.p1_GBM")

list.files()


library(dplyr)
library(data.table)
library(tidyverse)


data <- read.table("Supplementary Table 5. Bulk RNA-seq tpm matrix from ICB cohort..csv", 
                    sep = ",", header = T, check.names = F)
data <- distinct(.data = data, Description, .keep_all = T)
data <- column_to_rownames(.data = data, var = "Description")

cibersort <-read.table("Supplementary Table 18. CIBERSORTx abundance estimates from ICB cohort..csv", 
                       sep = ",", header = T, check.names = F)
cibersort <- cibersort[c(1:128),]

colnames(data) %in% cibersort$Mixture


library(dplyr)
library(tidyr)
library(tibble)

# 1. 从 data 矩阵中提取 IGSF3 的表达量，并转置为行为样本、列为基因的格式
# （注意：如果你的 data 还没做过 TPM/CPM/FPKM 归一化，建议先做 log2(TPM+1)）
if("IGSF3" %in% rownames(data)){
  igsf3_expr <- as.data.frame(t(data["IGSF3", , drop = FALSE])) %>% 
    rownames_to_column("Sample")
} else {
  stop("在 data 矩阵里没找到 IGSF3，请检查是否全是大写或存在拼写差异！")
}

# 2. 规范化 cibersort 的样本名
# 确保 Mixture 这一列是字符型，并重命名为 Sample 方便合并
cibersort_clean <- cibersort %>% 
  rename(Sample = Mixture)

# 3. 取样本交集并合并数据 (用 Sample 列对齐)
merged_data <- inner_join(igsf3_expr, cibersort_clean, by = "Sample")

# 4. 提取需要算相关性的细胞类型列表（排除第一列 Sample 和第二列 IGSF3）
cell_types <- setdiff(colnames(merged_data), c("Sample", "IGSF3"))

# 5. 循环计算 IGSF3 与每种细胞的相关系数 (R) 和显著性 p 值
cor_results <- lapply(cell_types, function(cell){
  res <- cor.test(merged_data$IGSF3, merged_data[[cell]], method = "spearman")
  data.frame(
    Cell_Type = cell,
    Correlation = res$estimate,
    p_value = res$p.value
  )
}) %>% bind_rows()

# 6. 给显著性打上星星标签 (*), 方便写论文
cor_results <- cor_results %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ "ns"
  )) %>%
  arrange(desc(Correlation)) # 按相关系数降序排序

library(ggplot2)
library(forcats)

ggplot(cor_results, aes(x = fct_reorder(Cell_Type, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  # 经典红蓝冷暖色调：正相关用红，负相关用蓝
  scale_fill_gradient2(low = "#4169E1", mid = "white", high = "#DC143C", midpoint = 0, limit = c(-1,1)) +
  # 在柱子上方或下方标注星号
  geom_text(aes(label = Significance, 
                y = Correlation + ifelse(Correlation >= 0, 0.03, -0.05)), 
            vjust = 0.5, size = 4, fontface = "bold") +
  coord_flip() + # 横过来排版最美观
  theme_bw() +
  theme(
    axis.text = element_text(size = 11, color = "black", face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  labs(
    x = "Immune Cell Types (CIBERSORT)",
    y = "Spearman Correlation Coefficient (R)",
    title = "Correlation between IGSF3 Expression & Immune Infiltration",
    fill = "Correlation (R)"
  )

# ggsave("IGSF3_Cibersort_Correlation.pdf", width = 7, height = 5)









