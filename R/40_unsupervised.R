library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)

# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")

# Normalization log2(fpkm+1)
log2fpkm <- cbind(fpkm_all[,1:7], log(fpkm_all[,8:ncol(fpkm_all)] + 1, base = 2))

# Remove genes without variance between patients
variances <- apply(log2fpkm[, 8:ncol(log2fpkm)], 1, var)
sd <- apply(log2fpkm[, 8:ncol(log2fpkm)], 1, sd)
q1 <- quantile(variances, na.rm = T)["25%"]
log2fpkm <- log2fpkm[-which(is.na(variances) | variances <= q1 | sd == 0), ]


# apply combat for batch effect

# parametric adjustment
combat_edata1 = ComBat(dat=as.matrix(log2fpkm[,8:ncol(log2fpkm)]), batch=p_all$Batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
# Batches by PCA
pca_rna <- prcomp(t(log2fpkm[,8:ncol(log2fpkm)]))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 
pca_rna <- prcomp(t(combat_edata1))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 


# comparison between replicates, correlation R and pval

# 11817
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_44", y = 'R4163_YZ_6',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_44", y = 'R4163_YZ_6',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 11840
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_45", y = 'R4163_YZ_28',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_45", y = 'R4163_YZ_28',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 12889
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_55", y = 'R4163_YZ_10',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_55", y = 'R4163_YZ_10',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 12890
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_6", y = 'R4163_YZ_27',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_6", y = 'R4163_YZ_27',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 12929
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_20", y = 'R4163_YZ_11',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_20", y = 'R4163_YZ_11',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 13034
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_8", y = 'R4163_YZ_12',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_8", y = 'R4163_YZ_12',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 13155
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_9", y = 'R4163_YZ_14',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_9", y = 'R4163_YZ_14',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 15002
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_40", y = 'R4163_YZ_24',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_40", y = 'R4163_YZ_24',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 



# Get the average of duplicated genes
result <- WGCNA::collapseRows(combat_edata1,
                          rowGroup=log2fpkm$Feature_gene_name,
                          rowID=rownames(combat_edata1),
                          method="function",
                          methodFunction=colMeans)

data <- data.frame(result$datETcollapsed) #use this for undupervised analysis
variances <- apply(data, 1, var)
sdv <- apply(data, 1, sd)


top_unsup <- data.frame(data) %>%
  mutate(gene =rownames(.),
    sdv = sdv,
    variances = variances) %>%
  arrange(desc(variances)) %>%
  select(gene) %>%
  head(100) %>%
  pull()

filtered_data <- data.frame(data) %>%
  mutate(gene =rownames(.)) %>%
  filter(gene %in% top_unsup) %>%
  select(.,gene, everything())

rownames(filtered_data) <- filtered_data$gene

heatmap(as.matrix(filtered_data[,2:ncol(filtered_data)]), scale='row')


cormat <- Hmisc::rcorr(t(as.matrix(filtered_data[,2:ncol(filtered_data)])))

heatmap(cormat$r)

# try removing unsignificant p vals from matrix


