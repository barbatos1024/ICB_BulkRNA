
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/Braun_2020_RCC/")

library(tidyverse)
library(data.table)

df <- read.table('Normalized_RNA-seq.csv', header = T, sep = ",")
df <- distinct(.data = df, gene_name, .keep_all = T)
colnames(df)
rownames(df) <- df$gene_name
df <- df[,-1]

pdata <- fread("pdata.csv", header = T, sep = ",")
colnames(pdata)
pdata <- pdata[,c(7,24)]
pdata <- na.omit(pdata)

colnames(df) %in% pdata$RNA_ID
table(pdata$Benefit)

pdata <- dplyr::arrange(pdata, Benefit)
sampleid <- pdata$RNA_ID

data <- df[,c(sampleid)]

group <- c(rep("CB", 89), rep("ICB", 120), rep("NCB", 102))

library(limma)

data_normalized <- normalizeBetweenArrays(data)

boxplot(data_normalized, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)")  

boxplot(data, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)")     


###箱线图
data_normalized <- as.data.frame(data_normalized)
GENE = data_normalized["IGSF3",]
z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("CB", 89), rep("ICB", 120), rep("NCB", 102))
z$type <- factor(z$type, levels = c("CB", "ICB", "NCB"))

# z$type = c(rep("CB", 209), rep("NCB", 102))
# z$type <- factor(z$type, levels = c("CB", "NCB"))

colnames(z)

library(ggpubr)


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
    comparisons = list(c("CB", "ICB"),
                       c("CB", "NCB"),
                       c("ICB", "NCB")),
    label.y = c(34.3, 34.8, 33.8), method = "t.test"  ### 或用"t.test" 进行配对非参数检验
  ) + ggtitle("IGSF3")

ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)
