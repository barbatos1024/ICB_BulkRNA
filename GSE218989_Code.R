
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE218989_NSCLC/")
getwd()

library(data.table)
library(tidyverse)

tpm <- read.table("GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt", header = T, row.names = 1)
colnames(tpm)

metadata <- read.table("GSE218989_metadata.csv", header = T, sep = ",")
metadata <- metadata[c(1:497),]   # 1:R   0:NR

samples <- colnames(tpm)
metadata_filtered <- metadata %>% 
  filter(PatientID %in% samples) %>%
  arrange(desc(Responder))

table(metadata_filtered$Responder) # 1:R_171   0:NR_184
patientid <- metadata_filtered$PatientID

tpm_ordered <- tpm[,c(patientid)]


### 创建分组 ###
group <- c(rep("R",171), rep("NR",184))

library(limma)

logtpm <- log2(tpm_ordered + 1)
logtpm_normalized <- normalizeBetweenArrays(logtpm)

boxplot(logtpm_normalized, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)", ylim = c(-5,12))  

boxplot(logtpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 (TPM + 1) Boxplot", ylab = "Log2 (TPM + 1)", ylim = c(-5,12))     


###箱线图
logtpm_normalized <- as.data.frame(logtpm_normalized)
GENE = logtpm["IGSF3",]
z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("R",171), rep("NR",184))

z$type <- factor(z$type, levels = c("R", "NR"))
colnames(z)

library(ggpubr)

ggplot(z, aes(type, IGSF3)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.23)) +
  geom_line(colour = "black", linetype = "11") +
  xlab("") +
  ylab(paste("Gene Expression log2(TPM + 1)")) + ### ylab(paste("Gene Expression log2 (FPKM + 1)", "IGSF3"))
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

ggsave("ISGF3_PD-1.png", width = 3, height = 5, dpi = 300)

