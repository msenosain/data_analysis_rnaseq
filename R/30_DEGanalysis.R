library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)
library(sva)
# Check BiocManager::valid()

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

    } else {
        data_t <- data.frame(t(normalized_data))
        mid_q <- quantile(data_t[,gene])[["50%"]]
        data_t$level <- "low"
        data_t$level[data_t[,gene]  > mid_q] <- "high"

    }

    meta_data <- cbind(metadata, level=data_t$level)
    meta_data$Batch <- factor(meta_data$Batch)
    meta_data$level <- factor(meta_data$level)

    meta_data
}

# For SLC7A11
gene_x <- 'KDM5D'
meta_data <- label_by_gene(vsd_mat, gene=gene_x, p_all, extremes_only=T)
meta_data <- meta_data %>% filter(level != "normal")
counts_all <- as.data.frame(counts_all[,meta_data$Vantage_ID])
pData_rnaseq <- pData_rnaseq %>% filter(pt_ID %in% meta_data$pt_ID)
vsd_mat <- vsd_mat[,meta_data$Vantage_ID]

dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = meta_data,
  design = ~Batch + level, tidy = F  
)

dds$level <- relevel(dds$level, ref = "low")
dds_deg <- DESeq(dds, parallel = F)

dds_level <- data.frame(results(dds_deg)) %>%
  mutate(gene=rownames(.))  %>%
  filter(padj < 0.001) %>%
  #dplyr::rename(gene = row) %>%
  arrange(desc(log2FoldChange))
  #arrange(padj)

top <- dds_level %>%
  head(100) %>%
  dplyr::select(gene) %>%
  pull()

#DH_gene_list <- read_excel("~/Downloads/DH_gene_list.xlsx")
#dh_genes <- DH_gene_list$Feature_gene_name

gene_list <- top
# Use normalized batch corrected collapsed data (vsd_mat) for plot
filtered_res <- data.frame(vsd_mat) %>%
  mutate(gene=rownames(.)) %>%
  filter(gene %in% gene_list) %>%
  #select(gene, all_of(low_high_patients)) %>%
  as.data.frame()
rownames(filtered_res) <- filtered_res$gene
filtered_res$gene <- NULL

# Plot

library(ComplexHeatmap)
ha = HeatmapAnnotation(

    CANARY = as.factor(pData_rnaseq$CANARY),
    KDM5D_level = as.factor(meta_data$level),

    simple_anno_size = unit(0.5, "cm")
)

    Heatmap(as.matrix(filtered_res), name = "mat", row_km = 2,
      heatmap_legend_param = list(color_bar = "continuous"), 
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8), top_annotation = ha, 
      column_split =as.factor(meta_data$level))
