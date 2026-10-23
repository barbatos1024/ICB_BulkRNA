
rm(list = ls())   # 无response信息

setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE195832_HNSCC/")

data <- read.table("GSE195832_TJU_featurecounts_symbol.tsv.gz", header = T, row.names = 1)
colnames(data)

pre <- c(seq(1,56,2))
pre
post = pre + 1
post

data <- data[,c(pre, post)]
colnames(data)

### 创建分组 ###
group <- c(rep("pre",28), rep("post", 28))

library(edgeR)

### 数据预处理
# 构建 DGEList 对象
dgelist <- DGEList(counts = data, group = group)
# 过滤 low count 数据 例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1) >= 2                ###要求至少在2个样本中表达量大于1
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]                  
# 标准化 以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

# 获取 log2 CPM 表达矩阵
log_cpm <- cpm(dgelist_norm, log = TRUE)

# 画箱型图
boxplot(log_cpm, las = 2, outline = FALSE, col = "lightblue",
        main = "Log2 CPM Boxplot", ylab = "Log2 CPM")        ########！！！！！！！


###箱线图
log_cpm_normalized <- as.data.frame(log_cpm)
GENE = log_cpm_normalized["IGSF3",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep("pre",28), rep("post", 28))

z$type <- factor(z$type, levels = c("pre", "post"))

z$pairinfo = c(rep(1:28, 2))

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
    comparisons = list(c("pre", "post")), paired = T,
    method = "t.test",   ### 或用"t.test" 进行配对非参数检验
    label.y = c(9)
  ) + ggtitle("IGSF3")

ggsave("ISGF3_PD-1.png", width = 3, height = 5, dpi = 300)
