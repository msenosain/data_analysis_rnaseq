source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')

# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/prcnts_by_pt.RData")

# Select only good quality rna samples
ls_preprocessed <- preprocess_rna(path_rnaseq = 'rnaseq.RData', correct_batch = T, correct_gender = T)
samples_id <- colnames(ls_preprocessed$vsd_mat)

tpm_all_cbs <- cbind("sample_ID"=tpm_all_cbs[,1], tpm_all_cbs[,samples_id])
tpm_all <- tpm_all[,samples_id]
fpkm_dcv_cbs <- cbind("sample_ID"=fpkm_dcv_cbs[,1], fpkm_dcv_cbs[,samples_id])
fpkm_dcv <- fpkm_dcv[,samples_id]


# write txt files for xcell and cibersort
write.table(tpm_all_cbs, file = "deconvolution/input/rna_only/rnaseq_tpm_cbs.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(tpm_all, file = "deconvolution/input/rna_only/rnaseq_tpm.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)

write.table(fpkm_dcv_cbs, file = "deconvolution/input/rna_only/rnaseq_fpkm_cbs.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(fpkm_dcv, file = "deconvolution/input/rna_only/rnaseq_fpkm.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)

# Deconvolution

# MCP counter
mcp_dcv <- MCPcounter::MCPcounter.estimate(tpm_all, featuresType = "HUGO_symbols")
write.table(mcp_dcv, file = "deconvolution/output/rna_only/MCP/mcp_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)
mcp_fpkm_dcv <- MCPcounter::MCPcounter.estimate(fpkm_dcv, featuresType = "HUGO_symbols")
write.table(mcp_fpkm_dcv, file = "deconvolution/output/rna_only/MCP/mcp_fpkm_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)


# quanTIseq
qts_dcv <- immunedeconv::deconvolute(tpm_all, "quantiseq")
write.table(qts_dcv, file = "deconvolution/output/rna_only/QTS/qts_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
qts_fpkm_dcv <- immunedeconv::deconvolute(fpkm_dcv, "quantiseq")
write.table(qts_fpkm_dcv, file = "deconvolution/output/rna_only/QTS/qts_fpkm_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

##################################

# match pt ids cytof and rna
pt_cytof <- samples_id[which(ls_preprocessed$p_all$pt_ID %in% pData_cytof$pt_ID)]
tpm_all_cbs <- cbind("sample_ID"=tpm_all_cbs[,1], tpm_all_cbs[,pt_cytof])
tpm_all <- tpm_all[,pt_cytof]
fpkm_dcv_cbs <- cbind("sample_ID"=fpkm_dcv_cbs[,1], fpkm_dcv_cbs[,pt_cytof])
fpkm_dcv <- fpkm_dcv[,pt_cytof]

# write txt files for xcell and cibersort
write.table(tpm_all_cbs, file = "deconvolution/input/rna_cytof/rnaseq_tpm_cbs.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(tpm_all, file = "deconvolution/input/rna_cytof/rnaseq_tpm.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)

write.table(fpkm_dcv_cbs, file = "deconvolution/input/rna_cytof/rnaseq_fpkm_cbs.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(fpkm_dcv, file = "deconvolution/input/rna_cytof/rnaseq_fpkm.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)

# Deconvolution

# MCP counter
mcp_dcv <- MCPcounter::MCPcounter.estimate(tpm_all, featuresType = "HUGO_symbols")
write.table(mcp_dcv, file = "deconvolution/output/rna_cytof/MCP/mcp_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)
mcp_fpkm_dcv <- MCPcounter::MCPcounter.estimate(fpkm_dcv, featuresType = "HUGO_symbols")
write.table(mcp_fpkm_dcv, file = "deconvolution/output/rna_cytof/MCP/mcp_fpkm_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)


# quanTIseq
qts_dcv <- immunedeconv::deconvolute(tpm_all, "quantiseq")
write.table(qts_dcv, file = "deconvolution/output/rna_cytof/QTS/qts_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
qts_fpkm_dcv <- immunedeconv::deconvolute(fpkm_dcv, "quantiseq")
write.table(qts_fpkm_dcv, file = "deconvolution/output/rna_cytof/QTS/qts_fpkm_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)


# CIBERSORT and xCell deconv done online