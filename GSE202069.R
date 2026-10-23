rm(list = ls())


getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE202069_HCC/")

library(data.table)
library(tidyverse)
pdata = fread("GSE202069_pdata.csv")
pdata <- na.omit(pdata)
pts <- pdata$Label

tpm <- fread("GSE202069_tpm.txt.gz", header = T)
tpm <- column_to_rownames(tpm, "V1")
colnames(tpm)

therapypts <- intersect(pts, colnames(tpm))
tpm <- tpm[,therapypts]

R <- paste0("R", 1:7)
NR <- paste0("NR", 1:9)
colnames(tpm) <- c(R, NR)

group_list = c(rep("R", 7), rep("NR", 9))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "R")

library(limma) 
data = normalizeBetweenArrays(tpm + 1)
#去除批次效应后
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-4,14))
#去除批次效应前
boxplot(tpm, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,16))

design = model.matrix(~ group_list)
colnames(design) <- levels(group_list) 
rownames(design) <- colnames(data) 

fit = lmFit(data, design)
fit = eBayes(fit)
colnames(fit)
allDiff = topTable(fit, coef = "NR", adjust = 'fdr', number = Inf) 
### allDiff = topTable(fit, coef = 2, adjust = 'fdr', number = Inf) 

allDiff$significant = "stable"
colnames(allDiff)
log2(1.5) #变化1.5倍
log2(2) #变化2倍
allDiff$significant[allDiff$logFC >= 0.585 & allDiff$P.Value < 0.05] = "up" #符合这些为up in NR
allDiff$significant[allDiff$logFC <= -0.585 & allDiff$P.Value < 0.05] = "down" #符合这些为down in NR

write.csv(allDiff, file = "diffgene.csv")


######## boxplot ##########
library(ggplot2)
library(ggpubr)
data <- as.data.frame(data)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('R', 7), rep("NR", 9))

z$type = factor(z$type, levels = c("R", "NR"))

colnames(z)
 

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2(TPM + 1)")) + 
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("R", "NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(6)
  ) + ggtitle("IGSF3")

ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)



