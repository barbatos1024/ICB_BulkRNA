rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE281729_HNSCC/")

library(data.table)
library(tidyverse)
library(dplyr)
library(ggpubr)

pdata <- as.data.frame(t(read.table("GSE281729_pdata.csv", header = T, sep = ",")))
colnames(pdata) <- pdata[1,]
pdata <- pdata[c(2:68),]

pdata$status <- paste(pdata$Treatment, pdata$Overall_Response, sep = "_")
table(pdata$status)

pdata$status <- factor(pdata$status, levels = c("Pre_MR", "Post_MR", "Pre_R", "Post_R",
                                                "Pre_CR", "Post_CR", "Pre_NR", "Post_NR"))
pdata <- pdata %>% arrange(status)

# Pre_MR:5
# Post_MR:4
# Pre_R:18
# Post_R:17
# Pre_CR:3
# Post_CR:3
# Pre_NR:10
# Post_NR:7

rownames(pdata)
pdata$status

### unmatched sample
unmatched <- c("TJ3_51_1", "TJ3_09_1", "TJ3_47_1", "TJ3_11_1", "TJ3_31_1", "TJ3_48_1", "TJ3_44_2")
pdata_filtered <- pdata %>% filter(!rownames(pdata) %in% unmatched)

count <- read.table("GSE281729_Mastrolonardo_etal_expressions_logTPM_Nivo-IDO_Bulk-RNAseq_Processed_File.txt", 
                    header = T, sep = "\t")
count <- column_to_rownames(.data = count, var = "X")
count <- count[,c(2:68)]

count <- count[,rownames(pdata_filtered)]
table(pdata_filtered$status)
colnames(count)

# Pre_MR:4
# Post_MR:4
# Pre_R:16
# Post_R:16
# Pre_CR:3
# Post_CR:3
# Pre_NR:7
# Post_NR:7

group_list = c(rep('Pre_MR', 4), rep("Post_MR", 4), rep("Pre_R", 16), rep("Post_R", 16), 
               rep('Pre_CR', 3), rep("Post_CR", 3), rep("Pre_NR", 7), rep("Post_NR", 7))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "Pre_MR")


library(limma) 
data = normalizeBetweenArrays(count)
#去除批次效应后
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,10))
#去除批次效应前data#去除批次效应前
boxplot(count, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,12))

###箱线图
data <- as.data.frame(data)
GENE = data["IGSF8",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('Pre_MR', 4), rep("Post_MR", 4), rep("Pre_R", 16), rep("Post_R", 16), 
           rep('Pre_CR', 3), rep("Post_CR", 3), rep("Pre_NR", 7), rep("Post_NR", 7))

z$type <- factor(z$type, levels = c("Pre_MR", "Post_MR", "Pre_R", "Post_R",
                                    'Pre_CR', "Post_CR", "Pre_NR", "Post_NR"))
table(pdata$Overall_Response)
colnames(z)


####### 匹配一下 ####### 4个组
data <- as.data.frame(data)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('Pre_MR', 4), rep("Post_MR", 4), rep("Pre_R", 16), rep("Post_R", 16), 
           rep('Pre_CR', 3), rep("Post_CR", 3), rep("Pre_NR", 7), rep("Post_NR", 7))

z$type <- factor(z$type, levels = c("Pre_MR", "Post_MR", 
                                    "Pre_R", "Post_R",
                                    'Pre_CR', "Post_CR",
                                    "Pre_NR", "Post_NR"))

z$pairinfo = c(rep(1:4, 2), rep(5:20, 2), rep(21:23, 2), rep(24:30, 2))

ggplot(z, aes(type, IGSF3, fill = type)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 TPM")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Pre_MR", "Post_MR"), c("Pre_R", "Post_R"),
                       c('Pre_CR', "Post_CR"), c("Pre_NR", "Post_NR")),
    paired = TRUE, method = "wilcox.test"   ### 或用"t.test" 进行配对非参数检验
  )

####### 匹配一下 ###### 2个组
data <- as.data.frame(data)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('Pre_R', 4), rep("Post_R", 4), rep("Pre_R", 16), rep("Post_R", 16), 
           rep('Pre_R', 3), rep("Post_R", 3), rep("Pre_NR", 7), rep("Post_NR", 7))

z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))

z$pairinfo = c(rep(1:4, 2), rep(5:20, 2), rep(21:23, 2), rep(24:30, 2))

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 TPM")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Pre_R", "Post_R"), c("Pre_NR", "Post_NR")),
    paired = TRUE, method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(6.5, 6.3)
  ) +
  stat_compare_means(
    comparisons = list(c("Post_R", "Post_NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7)
  ) + ggtitle("IGSF3")

ggsave("ISGF3_PD-1_IDO.png", width = 3, height = 5, dpi = 300)

###################### PD-1 mono ########################
colnames(pdata_filtered)
monopdata <- filter(.data = pdata_filtered, SecondDrug == "Placebo")

monocount <- count[,rownames(monopdata)]
table(monopdata$status)

group_list = c(rep('Pre_MR', 2), rep("Post_MR", 2), rep("Pre_R", 4), rep("Post_R", 4), 
               rep('Pre_CR', 1), rep("Post_CR", 1), rep("Pre_NR", 2), rep("Post_NR", 2))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "Pre_MR")


library(limma) 
monodata = normalizeBetweenArrays(monocount)
#去除批次效应后
boxplot(monodata, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,10))
#去除批次效应前data#去除批次效应前
boxplot(monocount, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-2,12))

####### 匹配一下 ####### 4个组
monodata <- as.data.frame(monodata)
GENE = monodata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('Pre_MR', 2), rep("Post_MR", 2), rep("Pre_R", 4), rep("Post_R", 4), 
           rep('Pre_CR', 1), rep("Post_CR", 1), rep("Pre_NR", 2), rep("Post_NR", 2))

z$type <- factor(z$type, levels = c("Pre_MR", "Post_MR", "Pre_R", "Post_R",
                                    'Pre_CR', "Post_CR", "Pre_NR", "Post_NR"))

z$pairinfo = c(rep(1:2, 2), rep(3:6, 2), rep(7, 2), rep(8:9, 2))

ggplot(z, aes(type, IGSF3, fill = type)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 TPM")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Pre_MR", "Post_MR"), c("Pre_R", "Post_R"),
                       c('Pre_CR', "Post_CR"), c("Pre_NR", "Post_NR")),
    paired = TRUE, method = "wilcox.test"   ### 或用"t.test" 进行配对非参数检验
  )

####### 匹配一下 ###### 2个组
monodata <- as.data.frame(monodata)
GENE = monodata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框

z$type = c(rep('Pre_R', 2), rep("Post_R", 2), rep("Pre_R", 4), rep("Post_R", 4), 
           rep('Pre_R', 1), rep("Post_R", 1), rep("Pre_NR", 2), rep("Post_NR", 2))


z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))

z$pairinfo = c(rep(1:2, 2), rep(3:6, 2), rep(7, 2), rep(8:9, 2))

ggplot(z, aes(type, IGSF3, fill = type)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 TPM")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("Pre_R", "Post_R"), c("Pre_NR", "Post_NR")),
    paired = TRUE, method = "wilcox.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(7, 7)
  )



