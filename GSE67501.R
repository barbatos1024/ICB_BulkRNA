###让报错变成英文 
Sys.setenv(LANGUAGE = "en")
#禁止chr转成factor
options(stringsAsFactors = FALSE)
###清空环境
rm(list = ls())

###设置路径
getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE67501_RCC/")

library(data.table)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library("illuminaHumanv4.db")

ls('package:illuminaHumanv4.db')

probe2symbol = toTable(illuminaHumanv4SYMBOL)

data <- fread("GSE67501_Non-normalized_data.txt.gz", header = T)
colnames(data)
data <- data[,c(1,2,4,6,8,10,12,14,16,18,20,22)]

matchid = probe2symbol[match(data$ID_REF, probe2symbol$probe_id), ] 
matchid = na.omit(matchid)

counts <- merge(x = data, y = matchid, by.x = names(data)[1], by.y = names(matchid)[1], all = F)
colnames(counts)
# R: 1 3 7 8
# NR: 2 4 5 6 9 10 11 12
counts <- counts[,c(13, 2:12)]

dataframe = distinct(counts, symbol, .keep_all = T) 
rownames(dataframe) = dataframe$symbol 
genesymbol <- dataframe$symbol 
dataframe <- as.data.frame(dataframe)
dataframe = dataframe[,c(2:12)] #去除多余的列
rownames(dataframe) = genesymbol

a1 = c(1,2,3,4,5,6,7)
a2 = c(8,9,10,11)

dataframe1 = dataframe[,c(a2,a1)]

group_list = c(rep('R',4),rep('NR',7))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "R")

library(limma) 

dataframe2 = log2(dataframe1)
dataframe3 = normalizeBetweenArrays(dataframe2)

boxplot(dataframe3, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(6, 16))
boxplot(dataframe2, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(5, 16))

design = model.matrix(~ group_list)
colnames(design) <- levels(group_list) #将差异分析矩阵的列名换成HC OA RA
rownames(design) <- colnames(dataframe3) #将差异分析矩阵的行名换成样本名


fit = lmFit(dataframe3, design)
fit = eBayes(fit) 

DEGs = topTable(fit, coef = 2, adjust = 'fdr', number = Inf) 

library(ggplot2)                         
library(ggrepel)                               

DEGs$significant = "stable" #添加新的significant一列，命名为stable
colnames(DEGs)
log2(2) #变化2倍
DEGs$significant[DEGs$logFC >= 0.585 & DEGs$P.Value < 0.05] = "up" #符合这些为up
DEGs$significant[DEGs$logFC <= -0.585 & DEGs$P.Value < 0.05] = "down" #符合这些为down

###箱线图
alldata3 <- as.data.frame(dataframe3)
GENE = alldata3["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("R", 4), rep("NR", 7))

z$type <- factor(z$type, levels = c("R", "NR"))

colnames(z)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Relative Gene Expression")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("R", "NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(11.3)
  ) + ggtitle("IGSF3")

ggsave("ISGF3_PD-1.png", width = 3, height = 5, dpi = 300)


