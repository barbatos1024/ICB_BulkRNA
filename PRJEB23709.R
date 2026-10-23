
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/PRJEB23709_SKCM/")


color1 = c("#53049d", "#7500aa", "#b17dd5", "#c587c7", 
           "#a0259f", "#d16ba5", "#dea7bc", "#f7b3ae")

color2 = c("#0f0a80", "#4b1475", "#6400a6", "#870fa1",
           "#aa2494", "#da596c", "#eb7954", "#f0bb3f")

color3 = c("#573299", "#6e3fc1", "#8a63cd", "#a78cd9", "#c4b3e8")

color4 = c("#1d136d", "#421b7b", "#600e86", "#9d2a7b",
           "#c44a60", "#ca5662", "#e5a544", "#fabdc0")

color5 = c("#360c6b", "#59187a", "#a3307a", "#d4466c", "#f46a57")

library(data.table)
library(tidyverse)
library(ggpubr)
library(org.Hs.eg.db)
library(AnnotationDbi)

count <- read.table("PD1-IPIPD1_counts.txt", header = T, sep = "\t")
colnames(count)

ensembl_ids <- count$ID

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

count$genesymbol <- gene_symbols
count <- distinct(.data = count, genesymbol, .keep_all = T)
count <- na.omit(count)
rownames(count) <- count$genesymbol
colnames(count)
count <- count[,-c(1,93)]

PD1 <- count[,c(41:90)]
PD1_CTLA4 <- count[,c(1:40,91)]

# write.csv(PD1, file = "PD1.csv")
# write.csv(PD1_CTLA4, file = "PD1_CTLA4.csv")


### PD1 ###
colnames(PD1)
pdata_pd1 <- read.table("Table S1. Clinical Summary of Patients Treated with Anti-PD-1 Monotherapy, Related to Table 1.csv", 
                        header = T, sep = ",")

# EDT_PR: 1,2,5,6 ___ R
# EDT_PD: 3,4,7,8 ___ NR
# EDT_SD: 9 ___ NR

# PRE_PR: 10,12,13,14,15,27,28,29,30,31,32,34,38,40,49   ___ R
# PRE_CR: 17,19,22,33                                     ___ R
# PRE_SD: 15,21,39,43,47,48                               ___ NR
# PRE_PD: 11,18,20,23,24,25,26,35,36,37,41,42,44,45,46,50   ___ NR

PD1_order <- PD1[,c(17,19,22,33, 10,12,13,14,15,27,28,29,30,31,32,34,38,40,49, 
                    1,2,5,6,
                    15,21,39,43,47,48, 11,18,20,23,24,25,26,35,36,37,41,42,44,45,46,50, 
                    9, 3,4,7,8)]


group <- c(rep("Pre_R", 19), rep("Post_R", 4), rep("Pre_NR", 22), rep("Post_NR", 5))

library(edgeR)
dgelist <- DGEList(counts = PD1_order, group = group)
keep <- rowSums(cpm(dgelist) > 1) >= 2               
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]                  
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

log_cpm <- cpm(dgelist_norm, log = TRUE)

boxplot(log_cpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 CPM Boxplot", ylab = "Log2 CPM")   

### 箱线图 ###
data <- as.data.frame(log_cpm)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 19), rep("Post_R", 4), rep("Pre_NR", 22), rep("Post_NR", 5))

z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))


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
    comparisons = list(c("Pre_R","Post_R"), c("Pre_NR", "Post_NR"), c("Post_R", "Post_NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(8.6, 8.2, 9)
  ) + ggtitle("IGSF3")

# ggsave("IGSF3_PD1.png", width = 4, height = 5, dpi = 300)


### PD1_CTLA4 ###
colnames(PD1_CTLA4)

pdata_pd1ctla4 <- read.table("Table S2. Clinical Summary of Patients Treated with Combined Anti-PD-1 and Anti-CTLA-4 Immunotherapy, Related to Table 1.csv", 
                             header = T, sep = ",")

# EDT_PR: 4,41 ___ R
# EDT_CR: 2,3,5 ___ R
# EDT_PD: 6,7,8 ___ NR
# EDT_SD: 1 ___ NR

# PRE_PR: 9,17,20,24,28,32,34,35,36,37,38   ___ R
# PRE_CR: 10,12,16,19,23,26,27,30,31,33     ___ R
# PRE_SD: 15,22,25,29,39  ___ NR
# PRE_PD: 11,13,14,18,21,40   ___ NR

PD1_CTLA4_order <- PD1_CTLA4[,c(10,12,16,19,23,26,27,30,31,33, 9,17,20,24,28,32,34,35,36,37,38,
                                2,3,5, 4,41,
                                15,22,25,29,39, 11,13,14,18,21,40,
                                1, 6,7,8)]

group <- c(rep("Pre_R", 21), rep("Post_R", 5), rep("Pre_NR", 11), rep("Post_NR", 4))

library(edgeR)
dgelist <- DGEList(counts = PD1_CTLA4_order, group = group)
keep <- rowSums(cpm(dgelist) > 1) >= 2               
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]                  
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

log_cpm <- cpm(dgelist_norm, log = TRUE)

boxplot(log_cpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 CPM Boxplot", ylab = "Log2 CPM")   

### 箱线图 ###
data <- as.data.frame(log_cpm)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 21), rep("Post_R", 5), rep("Pre_NR", 11), rep("Post_NR", 4))

z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))


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
    comparisons = list(c("Pre_R","Post_R"), c("Pre_NR", "Post_NR"), c("Post_R", "Post_NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(8.8, 8.7, 9.4)
  ) + ggtitle("IGSF3")

# ggsave("IGSF3_PD1_CTLA4.png", width = 4, height = 5, dpi = 300)


### Total ###
colnames(PD1_order)
group_PD1 <- c(rep("Pre_R", 19), rep("Post_R", 4), rep("Pre_NR", 22), rep("Post_NR", 5))
colnames(PD1_CTLA4_order)
group_PD1_CTLA4 <- c(rep("Pre_R", 21), rep("Post_R", 5), rep("Pre_NR", 11), rep("Post_NR", 4))

combined <- cbind(PD1_order[,c(1:19)], PD1_CTLA4_order[,c(1:21)],
                  PD1_order[,c(20:23)], PD1_CTLA4_order[,c(22:26)],
                  PD1_order[,c(24:45)], PD1_CTLA4_order[,c(27:37)],
                  PD1_order[,c(46:50)], PD1_CTLA4_order[,c(38:41)])
colnames(combined)
group <- c(rep("Pre_R", 40), rep("Post_R", 9), rep("Pre_NR", 33), rep("Post_NR", 9))

library(edgeR)
dgelist <- DGEList(counts = combined, group = group)
keep <- rowSums(cpm(dgelist) > 1) >= 2               
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]                  
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

log_cpm <- cpm(dgelist_norm, log = TRUE)

boxplot(log_cpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 CPM Boxplot", ylab = "Log2 CPM")  

### 箱线图 ###
data <- as.data.frame(log_cpm)
GENE = data["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 40), rep("Post_R", 9), rep("Pre_NR", 33), rep("Post_NR", 9))

z$type <- factor(z$type, levels = c("Pre_R", "Post_R", "Pre_NR", "Post_NR"))


ggplot(z, aes(type, IGSF3)) +
  geom_boxplot(outlier.shape = NA) +   # geom_point(size = 2, alpha = 0.5) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2 cpm")) +
  theme_bw() +
  ylim(0,12) +
  geom_jitter(aes(colour = group), shape = 16, position = position_jitter(0.2), size = 2) +
  scale_colour_manual(values = color5) +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(vjust = 1, hjust = 0.5, size = 14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),    # 去除网格线
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11)) +
  stat_compare_means(
    comparisons = list(c("Pre_R","Post_R"), c("Pre_NR", "Post_NR"), c("Post_R", "Post_NR")),
    method = "wilcox.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(8.8, 9.4, 10.6)
  ) 
# + ggtitle("PRJEB23709")

# ggsave("PRJEB23709_IGSF3_combined.svg")

