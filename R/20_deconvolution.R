source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')

# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
ls_preprocessed <- preprocess_rna(path_rnaseq = 'rnaseq.RData', correct_batch = T, correct_gender = T)

# write txt files for xcell and cibersort
write.table(tpm_all_cbs, file = "rnaseq_tpm_cbs.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(tpm_all, file = "rnaseq_tpm.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)

write.table(fpkm_dcv_cbs, file = "rnaseq_fpkm_cbs.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(fpkm_dcv, file = "rnaseq_fpkm.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)

# Deconvolution

# MCP counter
mcp_dcv <- MCPcounter::MCPcounter.estimate(tpm_all, featuresType = "HUGO_symbols")
write.table(mcp_dcv, file = "mcp_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)
mcp_fpkm_dcv <- MCPcounter::MCPcounter.estimate(fpkm_dcv, featuresType = "HUGO_symbols")
write.table(mcp_fpkm_dcv, file = "mcp_fpkm_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)


# quanTIseq
qts_dcv <- immunedeconv::deconvolute(tpm_all, "quantiseq")
write.table(qts_dcv, file = "qts_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
qts_fpkm_dcv <- immunedeconv::deconvolute(fpkm_dcv, "quantiseq")
write.table(qts_fpkm_dcv, file = "qts_fpkm_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)


# CIBERSORT and xCell deconv done online

# Read in everything
mcp_dcv <- read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/MCPcounter/mcp_dcv.txt", row.names=1)
qts_dcv <- read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/quanTIseq/qts_dcv.txt", row.names=1)
cbs_dcv <- t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/CIBERSORT/CIBERSORT.Output_Job8.txt", row.names=1))
xcell_dcv <- read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/xCell/xCell_rnaseq_tpm_xCell_0746040120.txt", row.names=1)
