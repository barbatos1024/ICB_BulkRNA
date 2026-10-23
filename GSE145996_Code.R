
getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE145996_SKCM/")

library(tidyverse)

fpkm <- read.table("GSE145996_Melanoma_Immune_FPKM.txt", header = T, sep = "\t")
fpkm <- distinct(.data = fpkm, Gene.Name, .keep_all = T)
rownames(fpkm) <- fpkm$Gene.Name
fpkm <- fpkm[,-c(1,2)]
colnames(fpkm)

# R: CR,PR 1-8
# NR: PD,SD 9-14

group <- c(rep("R", 8), rep("NR", 6))

library(limma)

fpkm_log <- log2(fpkm + 1)
fpkm_log_normalized <- normalizeBetweenArrays(fpkm_log)

boxplot(fpkm_log_normalized, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (FPKM + 1) Boxplot", ylab = "Log2 (FPKM + 1)", ylim = c(-1,11))  

boxplot(fpkm_log, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (FPKM + 1) Boxplot", ylab = "Log2 (FPKM + 1)", ylim = c(-2,12))     


###箱线图
fpkm_log_normalized <- as.data.frame(fpkm_log_normalized)
GENE = fpkm_log_normalized["IGSF3",]
z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("R", 8), rep("NR", 6))

z$type <- factor(z$type, levels = c("R", "NR"))
colnames(z)

library(ggpubr)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +   # geom_point(size = 2, alpha = 0.5) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 (FPKM + 1)")) +
  theme_classic() +
  geom_jitter(shape = 16, position = position_jitter(0.2), colour = "red") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14),
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("R", "NR")),
    paired = F, method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7.5)
    ) + ggtitle("IGSF3")


ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)

