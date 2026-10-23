
rm(list=ls())

library(GEOquery)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(ggplot2)
library(data.table)

getwd()
setwd("/Users/naturexv/Datasets//PD-1/BulkRNA/GSE91061/")

color1 = c("#53049d", "#7500aa", "#b17dd5", "#c587c7", 
           "#a0259f", "#d16ba5", "#dea7bc", "#f7b3ae")

color2 = c("#0f0a80", "#4b1475", "#6400a6", "#870fa1",
           "#aa2494", "#da596c", "#eb7954", "#f0bb3f")

color3 = c("#573299", "#6e3fc1", "#8a63cd", "#a78cd9", "#c4b3e8")

color4 = c("#1d136d", "#421b7b", "#600e86", "#9d2a7b",
           "#c44a60", "#ca5662", "#e5a544", "#fabdc0")

color5 = c("#360c6b", "#59187a", "#a3307a", "#d4466c", "#f46a57")

library(data.table)
expres = fread("GSE91061_BMS038109Sample.hg19KnownGene.fpkm.txt", data.table = F)
class(expres) #查看对象的数据结构
colnames(expres)

#将Entrez id转为gene symbol
library(org.Hs.eg.db)
##查看一下有哪些基因ID
keytypes(org.Hs.eg.db)
x1 = expres$V1
x1 = as.character(x1)
x2 = AnnotationDbi::select(org.Hs.eg.db, keys = x1, 
                           columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
data <- merge(x = x2, y = expres, by.x = 1, by.y = 1)
rownames(data) <- data$SYMBOL #有重复的基因名
data2 = distinct(data, SYMBOL, .keep_all = T) #根据gene_symbol，去除重复rownames；
rownames(data2) = data2$SYMBOL #依旧报错，能存在缺失值
data3 = na.omit(data2) #去除NA值
rownames(data3) = data3$SYMBOL #genesymbol作为行名
colnames(data3)
data4 = data3[,-c(1,2)] #去除多余的列

#########对临床样本信息pdata进行处理#########################
pdata = fread("pdata.txt", data.table = F, fill = T)
pdata <- pdata[c(3,5),]
rownames(pdata) <- pdata$Sample_title
pdata <- pdata[,-1]
pdata <- t(pdata)
pdata <- as.data.frame(pdata)
colnames(pdata) <- c("Treatment", "Response")

rownames(pdata)
### 提取配对样本的样本名 ###
pdatapair <- pdata[-c(7,8,9,14,15,16,21,32,33,40,47,48,49,52,53,58,59,74,85,88,97,98,109),]
length(rownames(pdatapair)) ### 去掉未配对的，共86个，43对

########### 提取存在配对的样本 ############ 
y = rownames(pdatapair)
y[1:10]
splitname = strsplit(y, split = "_", fixed = T) 

samplename1 = sapply(splitname, function(splitname){splitname[1]}) 
samplename2 = sapply(splitname, function(splitname){splitname[2]}) 
samplename = paste0(samplename1, "_", samplename2)

colnames(data4)[103] <- c("Pt92_Pre") # !!! pt92样本没有pre，自己修改

samplename <- intersect(samplename, colnames(data4))

data5 <- data4[,c(samplename)] ### 提取完成
colnames(data5)
###
pdatapair <- pdatapair[-c(5,6),]
# write.csv(pdatapair, file = "pdatapair.csv")
data6 <- data5[,-c(5,6)] # 去掉 UNK 两个样本

a1 <- c(1,9,15,18,27,32,36,40,46,47,52,63,66,70,71,73,77,81) #Pre_PD
a2 <- c(3,5,11,19,23,26,29,33,37,44,49,56,59,61,75) #Pre_SD
a3 <- c(7,13,22,42,53,58,68,79,84) #Pre_PRCR
a5 <- c(3,5,11,19,23,26,29,33,37,44,49,56,59,61,75,7,13,22,42,53,58,68,79,84) #Pre_CB 

b1 <- c(2,10,16,17,28,31,35,39,45,48,51,64,65,69,72,74,78,82) #On_PD, 18例
b2 <- c(4,6,12,20,24,25,30,34,38,43,50,55,60,62,76) #On_SD, 15例
b3 <- c(8,14,21,41,54,57,67,80,83) #On_PRCR, 9例
b5 <- c(4,6,12,20,24,25,30,34,38,43,50,55,60,62,76,8,14,21,41,54,57,67,80,83) #On_CB
length(b5)
### 治疗前 样本 提取 ###
Pre <- data6[,c(a1,a5)]
### 治疗后 样本 提取 ###
On <- data6[,c(b1,b5)]
### 配对的全部样本 ###
alldata <- data6[,c(a1,b1,a5,b5)] 
colnames(alldata)
# write.csv(alldata, file = "alldata_CB&PD.csv")

################################
alldata <- fread("alldata_CB&PD.csv")
alldata <- column_to_rownames(alldata, "V1")
alldata2 <- alldata
######################## 分组 ##########################
group_list = c(rep('Pre_PD', 18), rep("On_PD", 18), rep("Pre_CB", 24), rep("On_CB", 24))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "Pre_PD")

library(limma) 
alldata2 = log2(alldata + 1)
##除去批次效应
alldata3 = normalizeBetweenArrays(alldata2)
#去除批次效应后
boxplot(alldata3, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-1,9))
#去除批次效应前
boxplot(alldata2, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-1,10))

###箱线图
alldata3 <- as.data.frame(alldata3)
GENE = alldata3["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('Pre_PD', 18), rep("On_PD", 18), rep("Pre_CB", 24), rep("On_CB", 24))

z$type <- factor(z$type, levels = c("Pre_CB", "Pre_PD", "On_CB", "On_PD"))

colnames(z)


p2 = ggboxplot(z, x = "type", y = "IGSF3",
               width = 0.6, fill = "type",
               notch = F, palette = c("grey", "red", "blue", "yellow"),
               add = "jitter", shape = "type") + 
  labs(x = "", y = "Relative Expression", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 0, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(aes(group = type), label.x = 2, label.y = 8) +
  ggtitle("IGSF6") 

p2

####################################################
####### 配对 #######
alldata3 <- as.data.frame(alldata3)
GENE = alldata3["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('Pre_PD', 18), rep("On_PD", 18), rep("Pre_CB", 24), rep("On_CB", 24))

z$type <- factor(z$type, levels = c("Pre_CB", "On_CB", "Pre_PD", "On_PD"))

z$pairinfo = c(rep(1:18, 2), rep(19:42, 2))


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  geom_jitter(aes(colour = type), shape = 16, position = position_jitter(0.2), size = 2) +
  scale_colour_manual(values = color5) +
  xlab("") +
  ylab(paste("Gene Expression log2 (FPKM + 1)")) +
  theme_bw() +
  ylim(0, 10) +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("Pre_CB", "On_CB"), c("Pre_PD", "On_PD")), 
    paired = TRUE, method = "wilcox.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7.4, 7.7)
  ) +
  stat_compare_means(
    comparisons = list(c("On_CB", "On_PD")), 
    paired = F, method = "wilcox.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(8.9)
  ) 
# + ggtitle("GSE91061")

# ggsave("GSE91061_IGSF3_PD-1.svg")
