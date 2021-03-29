#---------------------------------------------------------------------------
# TF activity
#---------------------------------------------------------------------------

# Load functions and data
source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
environment_set()

ls_preprocessed <- preprocess_rna(path_rnaseq = '/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData', 
                                  correct_batch = T, lowvargenesrm = T, prot_coding_only = F)


# Function extracted from dorothea code
# https://github.com/saezlab/dorothea/blob/master/R/helpers.R#L17
dorothea2viper_regulons <- function(df) {
  regulon_list <- split(df, df$tf)
  viper_regulons <- lapply(regulon_list, function(regulon) {
    tfmode <- stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  return(viper_regulons)
}


result <- WGCNA::collapseRows(ls_preprocessed$vsd_mat,
                          rowGroup=ls_preprocessed$rna_all$Feature_gene_name,
                          rowID=rownames(ls_preprocessed$vsd_mat),
                          method="MaxMean")

data <- data.frame(result$datETcollapsed)

regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

regu <- dorothea2viper_regulons(regulons)
vpres_25 <- t(viper(data, regu, verbose = FALSE, minsize = 25))
rownames(vpres_25) <- p_all$pt_ID

vpres_4 <- t(viper(data, regu, verbose = FALSE, minsize = 4))
rownames(vpres_4) <- p_all$pt_ID


save(vsd_matTOP, vsd_matTOP_ENSEMBL, vsd_matTOP_clust, vsd_matTOP_clust_E, 
    clust_eigen, xcell_dcv, ed_dcv, vpres_25, vpres_4,
    file='/Users/senosam/Documents/Massion_lab/data_integration/RNA_data.Rdata')

