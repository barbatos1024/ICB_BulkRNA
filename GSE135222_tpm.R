rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE135222_NSCLC/")

library(data.table)
library(tidyverse)
library(dplyr)

tpm <- read.table("GSE135222_GEO_RNA-seq_omicslab_exp.tsv", header = T, sep = "\t")

library(biomaRt)
geneid <- tpm$gene_id
geneid_clean <- sub("\\..*", "", geneid)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = geneid_clean,
  mart = mart
)

tpm$ensembl_gene_id <- sub("\\..*", "", tpm$gene_id)  # 去掉版本号

data <- merge(gene_map, tpm, by = "ensembl_gene_id")
colnames(data)
data <- distinct(data, external_gene_name, .keep_all = T)
rownames(data) <- data$external_gene_name
data <- data[,-c(1,2,3)]


pdata <- read.table("GSE135222_pdata.csv", header = T, sep = ",")
pdata$cancer = "NSCLC"
pdata$sample = paste(pdata$cancer, pdata$PatientID, sep = "")

library(dplyr)

pdata_filtered <- pdata %>% 
  filter(sample %in% colnames(data)) %>%
  arrange(desc(benefit))

sampleid <- pdata_filtered$sample

data <- data[,c(sampleid)]   # R:8 NR:19


### 创建分组 ###
group <- c(rep("Pre_R", 8), rep("Pre_NR", 19)) 

library(limma)

logtpm <- log2(data + 1)
logtpm_normalized <- normalizeBetweenArrays(logtpm)

boxplot(logtpm_normalized, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)") 

boxplot(logtpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)")     


###箱线图
logtpm_normalized <- as.data.frame(logtpm_normalized)
GENE = logtpm_normalized["IGSF3",]
z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("Pre_R", 8), rep("Pre_NR", 19)) 

z$type <- factor(z$type, levels = c("Pre_R", "Pre_NR"))

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
    comparisons = list(c("Pre_R", "Pre_NR")),
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(5.5)
  ) + ggtitle("IGSF3")


ggsave("IGSF3.png", width = 3, height = 5, dpi = 300)
