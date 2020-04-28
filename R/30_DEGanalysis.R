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
p_all$Batch <- as.factor(p_all$Batch)
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = p_all,
  design = ~Batch, tidy = F  
)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)
vsd <- vst(dds)

# Remove batch effect
plotPCA(vsd, "Batch")
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Batch)
plotPCA(vsd, "Batch")
vsd_mat <- assay(vsd)

# Save pre-processed data
save(p_all, pData_rnaseq, counts_all, vsd_mat, file='dge_preprocessed.RData')

# Label samples
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/dge_preprocessed.RData")

label_by_gene <- function(normalized_data, gene, metadata, extremes_only=FALSE){
    if(extremes_only){
        data_t <- data.frame(t(normalized_data))
        first_q <- quantile(data_t[,gene])[["25%"]]
        third_q <- quantile(data_t[,gene])[["75%"]]
        data_t$level <- "normal"
        data_t$level[data_t[,gene] < first_q] <- "low"
        data_t$level[data_t[,gene] > third_q] <- "high"
        data_t$patient <- rownames(data_t)
        data_t$batch <- metadata$Batch

    } else {
        data_t <- data.frame(t(normalized_data))
        mid_q <- quantile(data_t[,gene])[["50%"]]
        data_t$level <- "low"
        data_t$level[data_t[,gene]  > mid_q] <- "high"
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

# For SLC7A11
meta_data <- label_by_gene(vsd_mat, gene='SLC7A11', p_all, extremes_only=F)

dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = meta_data,
  design = ~batch + level, tidy = F  
)

dds$level <- relevel(dds$level, ref = "low")
dds_deg <- DESeq(dds, parallel = F)
results(dds_deg)

dds_level <- data.frame(results(dds_deg)) %>%
  mutate(gene=rownames(.))  %>%
  #filter(padj < 0.001) %>%
  #dplyr::rename(gene = row) %>%
  arrange(desc(log2FoldChange))
  #arrange(padj)

top <- dds_level %>%
  head(200) %>%
  dplyr::select(gene) %>%
  pull()

DH_gene_list <- read_excel("~/Downloads/DH_gene_list.xlsx")
dh_genes <- DH_gene_list$Feature_gene_name


# Use normalized batch corrected collapsed data (vsd_mat) for plot
filtered_SLC7A11 <- data.frame(vsd_mat) %>%
  mutate(gene=rownames(.)) %>%
  filter(gene %in% top) %>%
  #select(gene, all_of(low_high_patients)) %>%
  as.data.frame()
rownames(filtered_SLC7A11) <- filtered_SLC7A11$gene
filtered_SLC7A11$gene <- NULL
colors <- c("red", "blue")
#levels <- as.factor(data_t %>% filter(level != "normal") %>% dplyr::select(level) %>% pull())
heatmap(as.matrix(filtered_SLC7A11), ColSideColors = colors[meta_data$level])






