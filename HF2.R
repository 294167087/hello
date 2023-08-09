if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")


library(tidyverse)
library(GEOquery)
library(limma)
library(tibble)
library(dplyr)


sample_info <- GSE116250_full_metadata <- read_csv("GSE116250_full_metadata.csv")
sample_info <- column_to_rownames(sample_info,var = "1")

gene_exp <- GSE116250_GeneLevel_Raw_data <- read_csv("GSE116250_GeneLevel_Raw_data.csv")
gene_exp <- column_to_rownames(gene_exp,var = "1")
gene_exp=log2(gene_exp+1)  

gene_exp=normalizeBetweenArrays(gene_exp)
boxplot(gene_exp,
        outline = FALSE,
        las=2)











gset <- getGEO("GSE21610",getGPL = FALSE)
gset <- gset$GSE21610_series_matrix.txt.gz
sample_info <- pData(gset)%>%
  dplyr::select(geo_accession,title)

write.table(sample_info,file = "D:/ZYT/R/HF/sample_info.txt",
            sep ="\t",
            quote = FALSE)



group_by(sample_info, group)%>%
  summarise(count= n())

gene_exp <- as.data.frame(exprs(gset))
gene_exp <- rownames_to_column(gene_exp,var = "probe")
gene_exp=log2(gene_exp+1)  

gene_exp=normalizeBetweenArrays(gene_exp)

boxplot(gene_exp,
        outline = FALSE,
        las=2)

save(gene_exp,sample_info,file = "D:/ZYT/R/HF/HF2.rdata")

##chayijiyinjianding limma/DEseq2/edgeR
##goujianshejijuzhen  design matrix
library(limma)
design <- model.matrix(~0 + sample_info$group)
colnames(design) <- levels(factor(sample_info$group))
rownames(design) <- rownames(sample_info)

##goujianbijiaojuzhen
contrasts <- makeContrasts(
  HN=HF-NF,
  levels = design)

##chayijiyin
fit <- lmFit(gene_exp,design)
fit <- contrasts.fit(fit,contrasts)
fit <- eBayes(fit)
de_result <- topTable(fit,
                      coef = "HN",
                      number = Inf)

gene_exp_tbl <- rownames_to_column(gene_exp,var = "probe")
gene_exp <- as.data.frame(gene_exp)

##tianjiajiyinxinxibiaohebiaodaliang tidyverse
tidyvese支持tibble,不是dataframe没有行名
de_result <- rownames_to_column(de_result,var = "probe")%>%
  mutate(direction = if_else(abs(logFC)<1 | P.Value>0.05,"ns",
                             if_else(logFC >= 1,"up","down")))
  
  de_result <- distinct(de_result,probe,.keep_all =TRUE)%>%
  left_join(gene_info,by = c("probe"= "ID"))%>%
  left_join(gene_exp_tbl,by= "probe")
  
temp <- select(de_result,-t,-B)
de_result <-  arrange(temp,desc(abs(logFC)))

write.table(de_result,file ="D:/ZYT/R/HF/HF2/de_result.txt",
            sep ="\t",
            quote = FALSE)

deg <- 
  #去除探针不对应基因的情况 或| 与& 非！
  filter(de_result,!is.na(gene_symbol))%>%
  #一个探针对应多个基因，保留第一个基因
  #方法一：一个基因对应多个探针，保留第一个
  distinct(gene_symbol,.keep_all =TRUE)
  #方法二：一个基因对应多个探针，保留abs(logFC)最大的
  group_by(ENTREZ_GENE_ID)%>%
  filter(abs(logFC) == max(abs(logFC)))%>%
  distinct(ENTREZ_GENE_ID,.keep_all = TRUE)

###提取差异基因列表
diffgene <- filter(deg,direction != "ns")

write.table(diffgene,file = "D:/ZYT/R/HF/HF2/diffgene.txt",
            sep ="\t",
            quote = FALSE)






library("ggplot2")
library("ggrepel")
hs_data <- read.delim("clipboard")
head(de_resultIS2)

ggplot(data = de_result, aes(x = logFC, y = -log10(P.Value)))+geom_point() 
de_result$threshold = as.factor(ifelse(de_result$P.Value < 0.05 & abs(de_result$logFC) >= 1.5, ifelse(de_result$logFC> 1.5 ,'Up','Down'),'NoSignifi'))
pvolcano <- ggplot(data = de_result, aes(x = logFC, y = -log10(P.Value), colour=threshold,label = probe)) +
  geom_point(alpha=0.4, size=3) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-5, 5)) +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)+
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential genes") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank())

###画热图
###提取所有IR2d,sham2d的样本编号
HF_NF_sample_info <- filter(sample_info,group%in% c("HF","NF"))

##提取用于画热图的差异表达基因的表格
temp<- slice(de_result,1:100)%>%
  select(Gene_Symbol,one_of(rownames(HF_NF_sample_info)))
temp <- as_tibble(temp)
temp <-  na.omit(temp)%>%
  distinct(Gene_Symbol,.keep_all =TRUE)
de_exp_top <-  column_to_rownames(temp,var = "Gene_Symbol")


cols <- list(group = c(HF= "red",
                       NF = "green"))



##热图
library(pheatmap)
IS2heatmap <- pheatmap(de_exp_top, 
                       annotation=select(HF_NF_sample_info,group), 
                       annotation_colors = cols,
                       color = colorRampPalette(c("green","white","red"))(50),
                       cluster_cols =F,
                       show_colnames = F,
                       show_rownames = F,
                       scale="row",
                       fontsize = 8,
                       fontsize_row=7,
                       fontsize_col=10)

