library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)
library(sva)
# Check BiocManager::valid()

# 0. Input: raw counts and fpkm
# 1. Pre-process:
#    - Remove duplicates (both), 
#    - Normalize and remove batch effect (fpkm)
# 2. Label samples: 
#    A. gene low/ high
#    B. clinical var
# 3. DEG




# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")

# Remove duplicates
dupl <- c('11817', '11840', '12889', '12890', '12929', '13034', '13155', '15002')
unique_vntg <- p_all %>%
  filter(., !(pt_ID %in% dupl &
    grepl("R4163",Vantage_ID))) %>%
  select(., Vantage_ID) %>%
  pull()

p_all <- p_all[which(p_all$Vantage_ID %in% unique_vntg),]
rna_all <- cbind(rna_all[,1:7], rna_all[,unique_vntg])
fpkm_all <- cbind(fpkm_all[,1:7], fpkm_all[,unique_vntg])
pData_rnaseq <- pData_rnaseq[!duplicated(pData_rnaseq$pt_ID),]
# pData_rnaseq$pt_ID == p_all$pt_ID # TRUE

# Remove low variance genes from counts
variances <- apply(rna_all[, 8:ncol(rna_all)], 1, var)
sdv <- apply(rna_all[, 8:ncol(rna_all)], 1, sd)
q1 <- quantile(variances, na.rm = T)["25%"]
idx <- which(is.na(variances) | variances <= q1 | sdv == 0)
rna_all <- rna_all[-idx, ]
result <- WGCNA::collapseRows(rna_all[, 8:ncol(rna_all)],
                          rowGroup=rna_all$Feature_gene_name,
                          rowID=rownames(rna_all),
                          method="function",
                          methodFunction=colMeans)
counts_all <- round(data.frame(result$datETcollapsed), digits = 0)

# Normalize
log2fpkm <- cbind(fpkm_all[,1:7], log(fpkm_all[,8:ncol(fpkm_all)] + 1, base = 2))
log2fpkm <- log2fpkm[-idx, ]
result <- WGCNA::collapseRows(log2fpkm[, 8:ncol(log2fpkm)],
                          rowGroup=log2fpkm$Feature_gene_name,
                          rowID=rownames(log2fpkm),
                          method="function",
                          methodFunction=colMeans)
log2fpkm <- data.frame(result$datETcollapsed)

# Remove batch effect
combat_log2fpkm = sva::ComBat(dat=as.matrix(log2fpkm), batch=p_all$Batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

## Batches by PCA
pca_rna <- prcomp(t(log2fpkm[,8:ncol(log2fpkm)]))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 
pca_rna <- prcomp(t(combat_log2fpkm))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 

# Save pre-processed data
save(p_all, pData_rnaseq, counts_all, combat_log2fpkm, file='dge_preprocessed.RData')

# Label samples
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/dge_preprocessed.RData")

label_by_gene <- function(normalized_data, gene, metadata, extremes_only=FALSE){
    if(extremes_only){
        data_t <- data.frame(t(normalized_data))
        mid_q <- quantile(data_t[,gene])[["50%"]]
        data_t$level <- "low"
        data_t$level[data_t[,gene]  > mid_q] <- "high"
        data_t$patient <- rownames(data_t)
        data_t$batch <- metadata$Batch

    } else {
        data_t <- data.frame(t(normalized_data))
        first_q <- quantile(data_t[,gene])[["25%"]]
        third_q <- quantile(data_t[,gene])[["75%"]]
        data_t$level <- "normal"
        data_t$level[data_t[,gene] < first_q] <- "low"
        data_t$level[data_t[,gene] > third_q] <- "high"
        data_t$patient <- rownames(data_t)
        data_t$batch <- metadata$Batch
    }

    meta_data <- data_t %>%
      select(patient, level, batch) %>%
      as.data.frame()
    meta_data$batch <- factor(meta_data$batch)
    meta_data$level <- factor(meta_data$level)

    meta_data
}


differential_expression_matrix <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = meta_data,
  design = ~batch + level, tidy = F  
)

differential_expression_matrix$level <- relevel(differential_expression_matrix$level, ref = "low")
differential_expression_analysis <- DESeq(differential_expression_matrix, parallel = F)
results(differential_expression_analysis)

differential_expression_result_level <- data.frame(results(differential_expression_analysis)) %>%
  mutate(gene=rownames(.))  %>%
  filter(padj < 0.001) %>%
  #dplyr::rename(gene = row) %>%
  arrange(desc(log2FoldChange))
  #arrange(padj)




# 1. Remove low variance genes from counts
variances <- apply(rna_all[, 8:ncol(rna_all)], 1, var)
sdv <- apply(rna_all[, 8:ncol(rna_all)], 1, sd)
q1 <- quantile(variances, na.rm = T)["25%"]
idx <- which(is.na(variances) | variances <= q1 | sdv == 0)
rna_all <- rna_all[-idx, ]
result <- WGCNA::collapseRows(rna_all[, 8:ncol(rna_all)],
                          rowGroup=rna_all$Feature_gene_name,
                          rowID=rownames(rna_all),
                          method="function",
                          methodFunction=colMeans)
counts_all <- round(data.frame(result$datETcollapsed), digits = 0)


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
result2 <- WGCNA::collapseRows(combat_edata1,
                          rowGroup=log2fpkm$Feature_gene_name,
                          rowID=rownames(combat_edata1),
                          method="function",
                          methodFunction=colMeans)

# 3. Label samples high normal and low for gene of interest
data <- data.frame(result2$datETcollapsed) 

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
count_data <- as.data.frame(counts_all[,low_high_patients])


# 4. DGE analysis

k <- which(meta_data$batch == 1)
meta_data$batch[k] = 'A'
meta_data$batch[-k] = 'B'
meta_data$batch <- factor(meta_data$batch)
meta_data$level <- factor(meta_data$level)
differential_expression_matrix <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = meta_data,
  design = ~batch + level, tidy = F  
)

differential_expression_matrix$level <- relevel(differential_expression_matrix$level, ref = "low")
differential_expression_analysis <- DESeq(differential_expression_matrix, parallel = F)
results(differential_expression_analysis)

differential_expression_result_level <- data.frame(results(differential_expression_analysis)) %>%
  mutate(gene=rownames(.))  %>%
  filter(padj < 0.001) %>%
  #dplyr::rename(gene = row) %>%
  arrange(desc(log2FoldChange))
  #arrange(padj)

differential_expression_result_level <- differential_expression_result_level[-c(1,2),] # removing SLC7A11 and its variant


top <- differential_expression_result_level %>%
  head(100) %>%
  dplyr::select(gene) %>%
  pull()

# Use normalized batch corrected collapsed data for plot

filtered_SLC7A11 <- data %>%
  mutate(gene=rownames(.)) %>%
  filter(gene %in% top) %>%
  select(gene, all_of(low_high_patients)) %>%
  as.data.frame()
rownames(filtered_SLC7A11) <- filtered_SLC7A11$gene
colors <- c("red", "blue")
levels <- as.factor(data_t %>% filter(level != "normal") %>% dplyr::select(level) %>% pull())
heatmap(as.matrix(filtered_SLC7A11[, 2:ncol(filtered_SLC7A11)]), ColSideColors = colors[levels])




