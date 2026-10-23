
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/OMIX013456_GC/")

list.files()

count <- read.table("OMIX013456-01_RNAseq_Counts.csv", header = T, sep = ",", row.names = 1)

MOS <- read.table("MOS_status.csv", header = T, sep = ",")
TCR <- read.table("TCR_Entropy.csv", header = T, sep = ",")
TCR$Sample_ID <- paste(TCR$Patient_ID, TCR$Time, sep = "_")
table(TCR$Sample_ID)

colnames(count) <- TCR$Sample_ID

library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)

genes <- rownames(count)
genes_clean <- sub("\\..*", "", genes)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = genes_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

count$gene_symbols <- gene_symbols
count <- distinct(.data = count, gene_symbols, .keep_all = T)
count <- na.omit(count)
rownames(count) <- count$gene_symbols
count <- count[,c(1:193)]
colnames(count)

### IGSF3 tcr correlation
ig_expr <- as.numeric(log2(count["IGSF3", ] + 1))
names(ig_expr) <- colnames(count)

df <- data.frame(
  Sample_ID = names(ig_expr),
  IGSF3 = ig_expr
)

df2 <- merge(df, TCR, by = "Sample_ID")

cor.test(df2$IGSF3, df2$TCR.Entropy, method = "pearson")

cor.test(df2$IGSF3, df2$TCR.Entropy, method = "spearman")

library(ggplot2)

ggplot(df2, aes(x = IGSF3, y = TCR.Entropy)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  theme_bw() +
  labs(
    x = "IGSF3 Expression",
    y = "TCR Entropy",
    title = "Correlation between IGSF3 and TCR Diversity"
  )

