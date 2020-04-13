library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)

# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")

# 1. Remove low variance genes from counts
variances <- apply(rna_all[, 8:ncol(rna_all)], 1, var)
sdv <- apply(rna_all[, 8:ncol(rna_all)], 1, sd)
q1 <- quantile(variances, na.rm = T)["25%"]
idx <- which(is.na(variances) | variances <= q1 | sdv == 0)
rna_all <- rna_all[-idx, ]

# 2. Use normalized counts batch corrected samples to find high normal and low expression of gene X
log2fpkm <- cbind(fpkm_all[,1:7], log(fpkm_all[,8:ncol(fpkm_all)] + 1, base = 2))
log2fpkm <- log2fpkm[-idx, ]

## parametric adjustment
combat_edata1 = ComBat(dat=as.matrix(log2fpkm[,8:ncol(log2fpkm)]), batch=p_all$Batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
## Batches by PCA
pca_rna <- prcomp(t(log2fpkm[,8:ncol(log2fpkm)]))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 
pca_rna <- prcomp(t(combat_edata1))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 

## Collapse duplicated genes
result <- WGCNA::collapseRows(combat_edata1,
                          rowGroup=log2fpkm$Feature_gene_name,
                          rowID=rownames(combat_edata1),
                          method="function",
                          methodFunction=colMeans)

# 3. Label samples high normal and low for gene of interest
data <- data.frame(result$datETcollapsed) 

## label low and high SLC7A11
data_t <- data.frame(t(data))
first_q <- quantile(data_t$SLC7A11)[["25%"]]
third_q <- quantile(data_t$SLC7A11)[["75%"]]
data_t$level <- "normal"
data_t$level[data_t$SLC7A11 < first_q] <- "low"
data_t$level[data_t$SLC7A11 > third_q] <- "high"
data_t$patient <- rownames(data_t)
data_t$batch <- p_all$Batch

meta_data <- data_t %>%
  select(patient, level,batch) %>%
  filter(level != "normal") %>%
  as.data.frame()
low_high_patients <- meta_data$patient
count_data <- as.data.frame(rna_all[,low_high_patients])\


# 4. DGE analysis

differential_expression_matrix <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = meta_data,
  design = ~batch + level, tidy = F  
)
differential_expression_matrix$level <- relevel(differential_expression_matrix$level, ref = "low")
differential_expression_analysis <- DESeq(differential_expression_matrix, parallel = F)
differential_expression_result_level <- as_tibble(results(differential_expression_analysis, tidy = T)) %>%
  filter(padj < 0.001) %>%
  dplyr::rename(gene = row) %>%
  # arrange(desc(log2FoldChange))
  arrange(padj)


