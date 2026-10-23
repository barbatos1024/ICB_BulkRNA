rm(list=ls())

### 加载r包
library(GEOquery)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(ggplot2)
##设置路径
getwd()
setwd("/Users/naturexv/Datasets/PD-1/BulkRNA/GSE78220")

####################################获取注释信息
library(data.table)
expres = fread("GSE78220_PatientFPKM.txt", data.table = F)  #获取注释信息 平台信息
class(expres) #查看对象的数据结构
colnames(expres)
rownames(expres) <- expres$Gene
expres <- expres[,-1]

pdata = fread("pdata.txt", data.table = F)
rownames(pdata) <- pdata$V1
pdata <- pdata[,-1]

###################################################查看分组信息
View(pdata)
colnames(expres)

a1 = c(1,6,9,10,12,14,16,17,18,19,23,24,25) #progressive disease
a2 = c(2,3,4,5,13,15,22,26,27,28,7,8,11,20,21) #response partical/completely

length(a2)
length(a1)
length(a3)
data = expres[,c(a2,a1)] #根据样本进行排序,response在前

respo <- paste0("respon", 1:15)
resis <- paste0("resis",1:13)

colnames(data) = c(respo, resis)

####################################################################
##这里我们创建分组信息，意思就是前四个是control，后四个是tumorexperienced
group_list = c(rep("respon", 15), rep("resis", 13))
##转化为因子 因子有顺序
group_list = factor(group_list) 
## 强制限定顺序
group_list <- relevel(group_list, ref = "respon")

library(limma) 
data2 <- log2(data + 0.001)
##除去批次效应
data3 = normalizeBetweenArrays(data2)
#去除批次效应后
boxplot(data3, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-15,14))

#去除批次效应前
boxplot(data2, outline = FALSE, notch = T, col = group_list, las = 2, ylim = c(-12,16))

##构建control和tumor差异分析的矩阵
design = model.matrix(~ group_list)
colnames(design) <- levels(group_list) #将差异分析矩阵的列名换成HC OA RA
rownames(design) <- colnames(data3) #将差异分析矩阵的行名换成样本名

##lmFit()：线性拟合模型构建
fit = lmFit(data3, design)
##eBayes()使用trend=TRUE对标准误差进行经验贝叶斯平滑，计算每个对比中每个基因的moderated t-statistic和log-odds。
fit = eBayes(fit) 
##topTable()给出一个最有可能在给定对比下差异表达的基因列表。
allDiff = topTable(fit, coef = 2, adjust = 'fdr', number = Inf)  #CR与PD的差异
?topTable
##topTable函数的coef参数，coef=2是指design的第2列，即tumour，即把tumour与normal进行对比
write.table(allDiff, file = "allDiff_ct.txt", sep = "\t", col.names = NA)
#col.names = NA 列名第一列空


#########################提取IGSF基因###########################################
###########################在gpl96中有的########################################
library(data.table)
IGSF = fread("IGSF.txt", header = T, sep = "\t")
colnames(IGSF)

IGSF = IGSF[,1]
GPL96 = splitID[,2]

GPL96ID = data.frame(GPL96) 
########################################
library(dplyr)
geneID = distinct(GPL96ID, GPL96, .keep_all = T) #根据gene_symbol，去除重复rownames；
?intersect
colnames(IGSF) = "genesymbol"
colnames(geneID) = "genesymbol"
IGSFgene = intersect(x = geneID, y = IGSF)
write.table(IGSFgene, file = "IGSFgenein96.txt", sep = "\t", col.names = NA)

library(ggvenn)
a1 <- as.matrix(IGSF)
b1 <- as.matrix(geneID)

a <- as.character(a1)
b <- as.character(b1)

y <- list("IGSF" = a, 
          "GPL" = b) 

p1 <- ggvenn(y, show_elements = F, fill_color = c("blue", "red"),
             label_sep = "\n", stroke_size = 1.5, set_name_size = 5,
             text_size = 3)
p1

IGSFgene = as.matrix(IGSFgene)
IGSFgene = as.character(IGSFgene)
###########################################################################################
##画个火山图
library(ggplot2)                         
library(ggrepel)                               

dataDIFF <- allDiff #这里计算tumor与control
dataDIFF$significant = "stable" #添加新的significant一列，命名为stable
colnames(dataDIFF)
log2(1.5) #变化1.5倍
log2(2) #变化2倍
dataDIFF$significant[dataDIFF$logFC >= 1 & dataDIFF$P.Value <0.05] = "up" #符合这些为up
dataDIFF$significant[dataDIFF$logFC <= -1 & dataDIFF$P.Value <0.05] = "down" #符合这些为down

write.csv(dataDIFF, file = "78220_2.csv")


p <- ggplot(dataCT,aes(logFC,-1*log10(adj.P.Val)))+xlim(-4,4)+ylim(0,3)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#EE7AE9","#f8961e"))+
  geom_hline(yintercept = 1.30,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-1, 1),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")

volcanoCT <- p + geom_label_repel(data=data_selected,
                                  aes(label=rownames(data_selected)))  
ggsave("volcanoCT.png", width = 8, height = 6.5, scale = 1)



#筛选基因
select.FPKM <- (data$AveExpr > 5) #平均表达量大于5
table(select.FPKM)

select.log2FC <- abs(data$logFC) > 0.585 #变化倍数大于1.5
table(select.log2FC)
select.qval <- (data$adj.P.Val< 0.05) #p<0.05
table(select.qval)

select.vec = (select.FPKM & select.log2FC & select.qval) #合并，同时满足
table(select.vec)

degs.list = as.character(rownames(data))[select.vec] 
#将差异基因赋值为degs.list并转换为character

label.deg = sample(degs.list, 20) #随机抽取20个
label.igsf = intersect(x = degs.list, y = IGSFgene) #抽取高表达中IGSF
write.table(label.igsf, file = "IGSF_diff.txt", sep = "\t", col.names = NA)

p=ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-5,5)+ylim(0,12)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#EE7AE9","#f8961e"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.585,0.585),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")

data_selected <- data[label.igsf, ] #根据label.deg提取data里的基因
ig = p + geom_label_repel(data=data_selected,
                          aes(label=rownames(data_selected)))
ig
ggsave("ig.png", width = 12, height = 10, scale = 1)



####################################################3
##画个热图
#创建数据框
label.kir <- c("KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", 
               "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DS1",
               "LILRB1", "LAG3", "PDCD1", "SIGLEC7", "CD300A", "CD96", "IL1RAPL1", "TIGIT",
               "HAVCR2", "KLRB1", "KLRC1")

KIRgene <- dataDIFF[label.kir,]
write.csv(KIRgene, "kir.csv", row.names = T)
##### kir原始表达 #####
kirgenes <- data[label.kir,]
write.csv(kirgenes, "kir_unlog.csv", row.names = T)

upgene <- filter(dataDIFF, significant == "up")
downgene <- filter(dataDIFF, significant == "down")
diffgene <- rbind(upgene, downgene)

data3 = as.data.frame(data3)  #转化为数据框
exprSet.map = data3[rownames(diffgene),]   #提取，提取20个基因/igsf基因

exprSet.map$number <- c(1:790)

mark <- c("KIR2DL3", "KIR3DL2")
mark
library(ComplexHeatmap)
#anno_mark()至少需要两个参数，其中at是原始数据矩阵的索引，标签是对应的文本。
lab = rowAnnotation(ano = anno_mark(at = c(332,457),
                                    labels = mark,
                                    labels_gp = gpar(fontsize = 12)))
library(colorRamp2)
col_fun = colorRamp2(c(-10, 2, 10), c("blue", "white", "red"))

Heatmap(exprSet.map,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 14),
        column_title = NULL,
        cluster_rows = T,
        cluster_columns = F,
        col = col_fun,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        right_annotation = lab,
        row_names_side = "left",
        name = "Exp")

############# 箱线图
data3 <- as.data.frame(data3)
GENE = data3["VSIG4",]

z = t(GENE) #转置
z = as.data.frame(z) #转化为数据框
z$type = c(rep('respo', 15), rep("resis", 13))

colnames(z)

p1 = ggboxplot(z, x = "type", y = "VSIG4",
               width = 0.6, fill = "type",
               notch = F, palette = c("grey", "red"),
               add = "jitter", shape = "type")

p1
###加个p值
p2 = ggboxplot(z, x = "type", y = "VSIG4",
               width = 0.6, fill = "type",
               notch = F, palette = c("grey", "red"),
               add = "jitter", shape = "type") + 
  labs(x = "", y = "Relative Expression", color = "") +
  theme(axis.title.x = element_text(vjust = 0, hjust = 0, size = 13), #调整x轴的title的角度，颜色，大小
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 0, color = "black", size = 11),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, angle = 0, color = "black", size = 11),
        legend.position = "none") +
  stat_compare_means(aes(group = type), label.x = 1.35, label.y = 8) +
  ggtitle("VSIG4") 

p2
ggsave("VSIG4.png", plot = p2, width = 3, height = 5, scale = 1)
####################################################
###箱线图
ggplot(z,aes(type,IGSF6,fill=type)) +
  geom_boxplot()

#######匹配一下
z$pairinfo=pairinfo=rep(1:10,3)

ggplot(z, aes(type,CTSC,fill=type)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","CTSC"))+
  theme_classic()+
  theme(legend.position = "none")

##########################################################################################
#####小提琴

ggplot(z,aes(type,IGSF6,fill=type)) +
  geom_violin()+ geom_jitter(shape=16, position=position_jitter(0.2))

library(ggsignif)
compaired <- list(c("HC", "OA", "RA"))

library(ggsci)
colnames(z)
igsf6_2=ggplot(z,aes(type,IGSF6,fill=type)) +
  geom_violin()+
  geom_signif(comparisons = compaired ,step_increase = 0.1,map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),test =t.test,size=2,textsize = 6)+
  geom_jitter(shape=16, position=position_jitter(0.2))

help(package="ggsci")

ggsave("igsf6_2.png", width = 6, height = 5, scale = 1)



#############################免疫浸润################################
data4 = normalizeBetweenArrays(data) #不需经log2处理,但需normalize
write.table(data4, "data4_imminfiltrated.txt", sep = "\t", col.names = NA)
write.table(dataframe7, "dataframe7.txt", sep = "\t", col.names = NA) #用log2处理的dataframe7会报错
## BiocManager::install("tidybulk")
## devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
source("Cibersort.R")
ciber <- CIBERSORT(sig_matrix = "LM22.txt",
                   mixture_file = "data4_imminfiltrated.txt",
                   perm = 1000,  #最好用1000，练习可用100
                   QN = TRUE)

ciber_data <- as.data.frame(ciber)
ciber_data$`NK cells activated`
ciber_data$`NK cells resting`
#######免疫浸润画图###########
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

colnames(ciber_data)
res <- data.frame(ciber_data[,1:22])%>%
  mutate(group = c(rep('R', 15), rep("NR", 13))) %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')

################################# Boxplot ##################################
res$group <- factor(res$group, levels = c("R", "NR"))

ggplot(res, aes(cell.type, value,fill = group)) + 
  geom_boxplot(outlier.shape = 21, color = "black") + 
  theme_bw() + ylim(0, 0.7) +
  labs(x = "", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 13, colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group), label = "p.format", size = 3, method = "kruskal.test") #这行为p值添加

immuneinfiltrate = ggplot(res, aes(cell.type, value,fill = group)) + 
  geom_boxplot(outlier.shape = 21, color = "black") + 
  theme_bw() + ylim(0, 0.7) +
  labs(x = "", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13, face = "italic", colour = 'black'),
        axis.text.y = element_text(face = "italic", size = 13, colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group), label = "p.format", size = 3, method = "kruskal.test") #这行为p值添加

ggsave("immuneinfiltrate.png", width = 20, height = 10, scale = 1)
############################################################################


############### rankplot ###############
### rank Fc ###
library(ggrepel)
colnames(KIRgene)
?arrange
KIRgene <- arrange(.data = KIRgene, desc(logFC)) # desc()降序
KIRgene$rankFC <- c(1:26)
rownames(KIRgene)

p1 <- ggplot(KIRgene, aes(x = rankFC, y = logFC)) + 
  geom_point(data = KIRgene[abs(KIRgene$logFC) > 0.585, ], size = 1, color = '#DC050C') +
  geom_point(data = KIRgene[abs(KIRgene$logFC) < 0.585,], 
             stroke = 0.5, size = 1, shape = 16, color = "grey", alpha = 1) +
  
  theme_classic() +
  geom_hline(yintercept = 0.585, linetype = 4, size = 0.3) +
  geom_hline(yintercept = -0.585, linetype = 4, size = 0.3) +
  theme(axis.title = element_text(colour = 'black', size = 15),
        axis.text = element_text(colour = 'black', size = 12)) +
  labs(x = 'Gene Rank', y = 'logFC')
p1

p2 <- p1 + geom_text_repel(data = KIRgene, direction = "both",
                           aes(x = rankFC, y = logFC,
                               label = rownames(KIRgene)))
p2
ggsave("waterfall_logFC.png", plot = p2, width = 28, height = 18, scale = 0.3)


p3 = p1 + geom_label_repel(data = KIRgene, force_pull = 5,
                           aes(x = rankFC, y = logFC, 
                               label = rownames(KIRgene)),
                           label.size = 0.15) # 标签框粗细
p3


### rank p.value ###
colnames(KIRgene)
?arrange
KIRgene <- arrange(.data = KIRgene, P.Value) # desc()降序
KIRgene$rankPV <- c(1:26)
rownames(KIRgene)

p1 <- ggplot(KIRgene, aes(x = rankPV, y = P.Value)) + 
  geom_point(data = KIRgene[KIRgene$P.Value < 0.05, ], size = 1, color = '#DC050C') +
  geom_point(data = KIRgene[KIRgene$P.Value > 0.05,], 
             stroke = 0.5, size = 1, shape = 16, color = "grey", alpha = 1) +
  
  theme_classic() +
  geom_hline(yintercept = 0.05, linetype = 4, size = 0.3) +
  theme(axis.title = element_text(colour = 'black', size = 15),
        axis.text = element_text(colour = 'black', size = 12)) +
  labs(x = 'Gene Rank', y = 'P.Value')
p1

p2 <- p1 + geom_text_repel(data = KIRgene, direction = "both",
                           aes(x = rankPV, y = P.Value,
                               label = rownames(KIRgene)))
p2
ggsave("waterfall_P.value.png", plot = p2, width = 28, height = 18, scale = 0.3)


p3 = p1 + geom_label_repel(data = kir, force_pull = 5,
                           aes(x = rankFC, y = logFC, 
                               label = rownames(kir)),
                           label.size = 0.15) # 标签框粗细
p3


