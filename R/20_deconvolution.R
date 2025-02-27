source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')

# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/prcnts_by_pt.RData")

main_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary'

# Select only good quality rna samples
ls_preprocessed <- preprocess_rna(path_rnaseq = file.path(main_path,'rnaseq.RData'), 
      correct_batch = T, correct_gender = T)
samples_id <- colnames(ls_preprocessed$vsd_mat)

tpm_all_cbs <- cbind("sample_ID"=tpm_all_cbs[,1], tpm_all_cbs[,samples_id])
tpm_all <- tpm_all[,samples_id]
fpkm_dcv_cbs <- cbind("sample_ID"=fpkm_dcv_cbs[,1], fpkm_dcv_cbs[,samples_id])
fpkm_dcv <- fpkm_dcv[,samples_id]

tpm_all2 <- tpm_all
rownames(tpm_all2) <- rna_all$Feature

result <- WGCNA::collapseRows(tpm_all2,
                          rowGroup=rna_all$Feature_gene_name,
                          rowID=rownames(tpm_all2),
                          method="MaxMean")

tpm_all_ED <- data.frame(result$datETcollapsed)


# write txt files for xcell and cibersort
write.table(tpm_all_cbs, file = file.path(main_path,"deconvolution/input/rna_only/rnaseq_tpm_cbs.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(tpm_all, file = file.path(main_path,"deconvolution/input/rna_only/rnaseq_tpm.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)
write.table(tpm_all_ED, file = file.path(main_path,"deconvolution/input/rna_only/rnaseq_tpm_ed.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)

write.table(fpkm_dcv_cbs, file = file.path(main_path,"deconvolution/input/rna_only/rnaseq_fpkm_cbs.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(fpkm_dcv, file = file.path(main_path,"deconvolution/input/rna_only/rnaseq_fpkm.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)

# Deconvolution

# MCP counter
mcp_dcv <- MCPcounter::MCPcounter.estimate(tpm_all, featuresType = "HUGO_symbols")
write.table(mcp_dcv, file = file.path(main_path,"deconvolution/output/rna_only/MCP/mcp_dcv.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)
mcp_fpkm_dcv <- MCPcounter::MCPcounter.estimate(fpkm_dcv, featuresType = "HUGO_symbols")
write.table(mcp_fpkm_dcv, file = file.path(main_path,"deconvolution/output/rna_only/MCP/mcp_fpkm_dcv.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)


# quanTIseq
qts_dcv <- immunedeconv::deconvolute(tpm_all, "quantiseq")
write.table(qts_dcv, file = file.path(main_path,"deconvolution/output/rna_only/QTS/qts_dcv.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
qts_fpkm_dcv <- immunedeconv::deconvolute(fpkm_dcv, "quantiseq")
write.table(qts_fpkm_dcv, file = file.path(main_path,"deconvolution/output/rna_only/QTS/qts_fpkm_dcv.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)

##################################

# match pt ids cytof and rna
pt_cytof <- samples_id[which(ls_preprocessed$p_all$pt_ID %in% pData_cytof$pt_ID)]
tpm_all_cbs <- cbind("sample_ID"=tpm_all_cbs[,1], tpm_all_cbs[,pt_cytof])
tpm_all <- tpm_all[,pt_cytof]
fpkm_dcv_cbs <- cbind("sample_ID"=fpkm_dcv_cbs[,1], fpkm_dcv_cbs[,pt_cytof])
fpkm_dcv <- fpkm_dcv[,pt_cytof]

# write txt files for xcell and cibersort
write.table(tpm_all_cbs, file = file.path(main_path,"deconvolution/input/rna_cytof/rnaseq_tpm_cbs.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(tpm_all, file = file.path(main_path,"deconvolution/input/rna_cytof/rnaseq_tpm.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)

write.table(fpkm_dcv_cbs, file = file.path(main_path,"deconvolution/input/rna_cytof/rnaseq_fpkm_cbs.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(fpkm_dcv, file = file.path(main_path,"deconvolution/input/rna_cytof/rnaseq_fpkm.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)

# Deconvolution

# MCP counter
mcp_dcv <- MCPcounter::MCPcounter.estimate(tpm_all, featuresType = "HUGO_symbols")
write.table(mcp_dcv, file = file.path(main_path,"deconvolution/output/rna_cytof/MCP/mcp_dcv.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)
mcp_fpkm_dcv <- MCPcounter::MCPcounter.estimate(fpkm_dcv, featuresType = "HUGO_symbols")
write.table(mcp_fpkm_dcv, file = file.path(main_path,"deconvolution/output/rna_cytof/MCP/mcp_fpkm_dcv.txt"), sep = "\t",
            row.names = TRUE, quote = FALSE)


# quanTIseq
qts_dcv <- immunedeconv::deconvolute(tpm_all, "quantiseq")
write.table(qts_dcv, file = file.path(main_path,"deconvolution/output/rna_cytof/QTS/qts_dcv.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
qts_fpkm_dcv <- immunedeconv::deconvolute(fpkm_dcv, "quantiseq")
write.table(qts_fpkm_dcv, file = file.path(main_path,"deconvolution/output/rna_cytof/QTS/qts_fpkm_dcv.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)


# EpiDISH with Vera's signature
library(EpiDISH)
signature <- as.matrix(read.table(file.path(main_path,'deconvolution/gene150_32.txt'), 
                                  header = TRUE, row.names = 1))
t1 <- read.table(file.path(main_path,'deconvolution/input/rna_only/rnaseq_tpm_ed.txt'), 
                 header = TRUE, row.names = 1)

res <- epidish(t1, signature, method = "RPC")
Fres <- as.data.frame(res$estF)
Fres <- cbind(Sample = rownames(Fres), Fres)
Fres[, -1] <- round(Fres[, -1], 3)

write.table(Fres, file.path(main_path, 'deconvolution/output/rna_only/epiDISH/ed_dcv.txt'), 
            sep = "\t", quote = F, row.names = F, col.names = T)


# CIBERSORT and xCell deconv done online