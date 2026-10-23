
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE244982")

library(data.table)
anno <- as.data.frame(t(read.table("GSE244982_series_matrix.txt", row.names = 1, header = T)))

exprs <- read.table("GSE244982_ProcessedData_bulkRNAseq.txt.gz", row.names = 1, header = T) #log2(fpkm + 1)
colnames(exprs)
anno$Sample_characteristics_ch0

a1 <- c(1:19, 23, 28) # CTLA-4_res
a2 <- c(20:22, 24:27, 29:41) # PD-1_res

exprs <- exprs[,c(a1,a2)] # 根据ctla4，pd-1耐药进行分组

group_list = c(rep('CTLA-4_res', 21), rep("PD-1_res", 20))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "CTLA-4_res")

library(limma)

##除去批次效应
data = normalizeBetweenArrays(exprs)
#去除批次效应后
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-1,9))
#去除批次效应前
boxplot(exprs, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-1,10))

###箱线图
library(ggpubr)

data <- as.data.frame(data)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('CTLA-4_res', 21), rep("PD-1_res", 20))

z$type <- factor(z$type, levels = c("CTLA-4_res", "PD-1_res"))

colnames(z)

p1 = ggboxplot(z, x = "type", y = "IGSF3",
               width = 0.6, fill = "type",
               notch = F, palette = c("grey", "red", "blue", "yellow"),
               add = "jitter", shape = "type") + 
  labs(x = "", y = "Gene Expression log2 (FPKM + 1)", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 0, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(aes(group = type), label.x = 1.3, label.y = 8) +
  ggtitle("IGSF3") 

p1

### 
ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +   # geom_point(size = 2, alpha = 0.5) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 (FPKM + 1)")) +
  theme_classic() +
  geom_jitter(shape = 16, position = position_jitter(0.2), colour = "red") +
  ylim(0, 8) +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("CTLA-4_res", "PD-1_res")), 
    paired = F, method = "wilcox.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7.5)
  )

ggsave("IGSF3_PD-1_CTLA4.png", width = 4, height = 5, dpi = 300)


