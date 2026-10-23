
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE294272_SKCM")

list.files()

library(edgeR)
library(ggpubr)
library(ggplot2)

counts <- read.table("GSE294272_counts.csv", header = T, sep = ",", row.names = 1, check.names = F)
metadata <- read.table("GSE294272_metadata.csv", header = T, sep = ",", check.names = F)

### untreated: 1,2,3,5,7,9,10,12,14,18,24,25,26,27             14
### response: 8,16,20,22,20                                    5
### resistance: 4,6,11,13,15,17,19,21,23,28                    10  
colnames(counts)

count <- counts[ , c(1,2,3,5,7,9,10,12,14,18,24,25,26,27, 8,16,20,22,20, 4,6,11,13,15,17,19,21,23,28)]
count_metadata <- data.frame(
  group = c(rep("Untrt", 14), rep("Response", 5), rep("Resistance", 10)),
  type = c(rep("LN", 29)),
  row.names = colnames(count)
)

all(rownames(count_metadata) %in% colnames(count))    # 检查是否一致
all(rownames(count_metadata) == colnames(count))    # 检查顺序一致

group = c(rep("Untrt", 14), rep("Response", 5), rep("Resistance", 10))

dgelist <- DGEList(counts = count, group = group)
keep <- rowSums(cpm(dgelist) > 1) >= 2               
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]                  
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

log_cpm <- cpm(dgelist_norm, log = TRUE)

data <- as.data.frame(log_cpm)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Untrt", 14), rep("Response", 5), rep("Resistance", 10))

z$type <- factor(z$type, levels = c("Untrt", "Response", "Resistance"))

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot(outlier.shape = NA) +   # geom_point(size = 2, alpha = 0.5) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 cpm")) +
  theme_classic() +
  geom_jitter(shape = 16, position = position_jitter(0.2), colour = "red") +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14),
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("Response","Resistance")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(10)
  ) + ggtitle("IGSF3")

# ggsave("IGSF3_PD1_CTLA4.png", width = 4, height = 5, dpi = 300)