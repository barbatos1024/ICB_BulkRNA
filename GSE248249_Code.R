rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/GSE248249_NSCLC/")

library(ggpubr)
library(data.table)

data <- read.table("GSE248249_nsclc_pd-1.csv", header = T, sep = ",", row.names = 1)

colnames(data)
data <- data[,c(1,7,10,12,16,19,25,27,33,35,37,40,42, # pre
                2,3,4,5,6,8,9,11,13,14,15,17,18,20,21,22,23,24,26,28,29,30,31,32,34,36,38,39,41)] #post

group_list <- c(rep("ICB_pre", 13), rep("ICB_post", 29))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "ICB_pre")

library(limma)
data_normalized <- normalizeBetweenArrays(data)

#去除批次效应后
boxplot(data_normalized, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(0, 10))
#去除批次效应前
boxplot(data, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(0, 10))

### 箱线图
data_normalized <- as.data.frame(data_normalized)
GENE = data_normalized["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("ICB_pre", 13), rep("ICB_post", 29))

z$type <- factor(z$type, levels = c("ICB_pre", "ICB_post"))

colnames(z)

### 加个p值

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
    comparisons = list(c("ICB_pre", "ICB_post")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(5.3)
  ) + ggtitle("IGSF3")

ggsave("ISGF3_PD-1.png", width = 2.5, height = 5, dpi = 300)


############# 配对样本 #############

data <- read.table("GSE248249_nsclc_pd-1.csv", header = T, sep = ",", row.names = 1)
data_paired <- data[,c(1,2,7,8,10,11,12,13,16,17,19,20,25,26,27,28,33,34,35,36,37,38,40,41,42,31)]
data_paired <- data_paired[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,   #pre
                              2,4,6,8,10,12,14,16,18,20,22,24,26)]   #post

group_list <- c(rep("ICB_pre", 13), rep("ICB_post", 13))
group_list = factor(group_list) 
group_list <- relevel(group_list, ref = "ICB_pre")

library(limma)
data_paired_normalized <- normalizeBetweenArrays(data_paired)

#去除批次效应后
boxplot(data_paired_normalized, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(0, 10))
#去除批次效应前
boxplot(data_paired, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(0, 10))

### 箱线图
data_paired_normalized <- as.data.frame(data_paired_normalized)
GENE = data_paired_normalized["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("ICB_pre", 13), rep("ICB_post", 13))

z$type <- factor(z$type, levels = c("ICB_pre", "ICB_post"))

colnames(z)
z$pairinfo = rep(1:13, 2)


###加个p值
library(ggpubr)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.5) +
  geom_line(aes(group = pairinfo), colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Relative Gene Expression")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
  theme_classic()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 0, hjust = 0.5, size = 14), #调整y轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(
    comparisons = list(c("ICB_pre", "ICB_post")),
    paired = TRUE, method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(4.5)
  ) + ggtitle("IGSF3")

ggsave("ISGF3_PD-1_paired.png", width = 2.5, height = 5, dpi = 300)



