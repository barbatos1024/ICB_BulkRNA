rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/Miya_PD-1_RCC_bulk/")

library(data.table)
library(tidyverse)
library(ggpubr)

normalizeddata <- read.table("PD-1_RCC_normalized.csv", header = T, sep = ",")
normalizeddata <- distinct(normalizeddata, Description, .keep_all = T)
rownames(normalizeddata) <- normalizeddata$Description
normalizeddata <- normalizeddata[,-c(1,2)]


pdata <- read.table("pdata.csv", header = T, sep = ",")
colnames(normalizeddata) %in% pdata$SUBJECT_ID 

samples <- colnames(normalizeddata)
pdata <- pdata %>% 
  filter(SUBJECT_ID %in% samples) %>%
  arrange(Response)
table(pdata$Response)

x = pdata$SUBJECT_ID

normalizeddata <- normalizeddata[,c(x)]
normalizeddata <- normalizeddata[,-c(89,90)]


group_list <- c(rep("CR", 9), rep("PD", 30), rep("PR", 22), rep("SD", 27))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "CR")


library(limma)

logdata = log2(normalizeddata + 1)

nordata = normalizeBetweenArrays(logdata)
#去除批次效应后
boxplot(nordata, outline = FALSE, notch = T, col = group_list, las = 2)
#去除批次效应前
boxplot(logdata, outline = FALSE, notch = T, col = group_list, las = 2)

### 箱线图
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("CR", 9), rep("PD", 30), rep("PR", 22), rep("SD", 27))

z$type <- factor(z$type, levels = c("PD", "SD", "PR", "CR"))

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
    comparisons = list(c('PD', 'SD'),
                       c('PD', 'PR'),
                       c('PD', 'CR')), 
    label.y = c(4.4, 4.7, 5.0), method = "t.test"
  ) + ggtitle("IGSF3")

ggsave("IGSF3.PNG", width = 4, height = 5, dpi = 300)

###############################
######################PD SD CR/PR
############################################
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("CRPR", 9), rep("PD", 30), rep("CRPR", 22), rep("SD", 27))

z$type <- factor(z$type, levels = c("PD", "SD", "CRPR"))

colnames(z)

p1 = ggboxplot(z, x = "type", y = "IGSF3",
               width = 0.6, fill = "type",
               notch = F,
               add = "jitter", shape = "type") + 
  labs(x = "", y = "Gene Expression log2(TPM + 1)", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 0, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c('PD', 'SD'),
                                        c('PD', 'CRPR')), 
                     label.y = c(4.4, 4.7), method = "t.test") +
  ggtitle("IGSF3") 

p1

###############################
######################NCB(PD) CB(SD CR/PR)
############################################
nordata <- as.data.frame(nordata)
GENE = nordata["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("CB", 9), rep("NCB", 30), rep("CB", 22), rep("CB", 27))

z$type <- factor(z$type, levels = c("NCB", "CB"))

colnames(z)

p1 = ggboxplot(z, x = "type", y = "IGSF3",
               width = 0.6, fill = "type",
               notch = F,
               add = "jitter", shape = "type") + 
  labs(x = "", y = "Gene Expression log2(TPM + 1)", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 0, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c('NCB', 'CB')), 
                     label.y = c(5), method = "t.test") +
  ggtitle("IGSF3") 

p1
