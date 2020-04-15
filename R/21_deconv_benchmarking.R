library(readr)
library(dplyr)
library(tidyr)

# Load deconvolution data
mcp_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/MCPcounter/mcp_dcv.txt", row.names=1)))
qts_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/quanTIseq/qts_dcv.txt", row.names=1)))
cbs_dcv <- read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/CIBERSORT/CIBERSORT.Output_Job8.txt", row.names=1)
xcell_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/xCell/xCell_rnaseq_tpm_xCell_0746040120.txt", row.names=1)))

# Load RNA seq data 
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")

# Load CyTOF data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/prcnts_by_pt.RData")

prcnt_pt <- prcnt_pt %>%
  tibble::rownames_to_column('pt_ID') %>%
  mutate(., Epithelial = rowSums(prcnt_pt[,9:ncol(prcnt_pt)])) %>%
  select(., Epithelial, everything()) %>%
  select(., -colnames(prcnt_pt)[8:ncol(prcnt_pt)]) %>%
  tibble::column_to_rownames('pt_ID')


# remove duplicates from rna
dupl <- c('11817', '11840', '12889', '12890', '12929', '13034', '13155', '15002')
unique_vntg <- p_all %>%
  filter(., !(pt_ID %in% dupl &
    grepl("R4163",Vantage_ID))) %>%
  select(., Vantage_ID) %>%
  pull()

unique_pt <- p_all %>%
  filter(., Vantage_ID %in% unique_vntg) %>%
  select(., pt_ID) %>%
  pull()

mcp_dcv <- mcp_dcv[unique_vntg,]
qts_dcv <- qts_dcv[unique_vntg,]
cbs_dcv <- cbs_dcv[unique_vntg,]
xcell_dcv <- xcell_dcv[unique_vntg,]

# row names by pt id
rownames(mcp_dcv) <- unique_pt
rownames(qts_dcv) <- unique_pt
rownames(cbs_dcv) <- unique_pt
rownames(xcell_dcv) <- unique_pt

# match pt ids cytof and rna
pt_cytof <- pData_cytof$pt_ID[which(pData_cytof$pt_ID %in% unique_pt)]
mcp_dcv <- mcp_dcv[pt_cytof,]
qts_dcv <- qts_dcv[pt_cytof,]
cbs_dcv <- cbs_dcv[pt_cytof,]
xcell_dcv <- xcell_dcv[pt_cytof,]
prcnt_pt <- prcnt_pt[pt_cytof,]

# Load excel cell types names
mcp_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "mcp")
qts_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "qts")
cbs_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "cbs")
xcell_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "xcell")

mcp_ctnames[]<- lapply(mcp_ctnames, function(x) gsub('\\+|\\-|\\s', '.', x))
qts_ctnames[]<- lapply(qts_ctnames, function(x) gsub('\\+|\\-|\\s', '.', x))
cbs_ctnames[]<- lapply(cbs_ctnames, function(x) gsub('\\+|\\-|\\s', '.', x))
xcell_ctnames[]<- lapply(xcell_ctnames, function(x) gsub('\\+|\\-|\\s', '.', x))

# match names
mcp_ctnames_v <- as.vector(as.matrix(mcp_ctnames))
mcp_ctnames_v <- mcp_ctnames_v[!is.na(mcp_ctnames_v)]

mcp_ctnames_v <- as.vector(as.matrix(mcp_ctnames))
mcp_ctnames_v <- mcp_ctnames_v[!is.na(mcp_ctnames_v)]

mcp_ctnames_v <- as.vector(as.matrix(mcp_ctnames))
mcp_ctnames_v <- mcp_ctnames_v[!is.na(mcp_ctnames_v)]

mcp_ctnames_v <- as.vector(as.matrix(mcp_ctnames))
mcp_ctnames_v <- mcp_ctnames_v[!is.na(mcp_ctnames_v)]


somefunction <- function(ctnames_df, deconv_results){
    ctnames_v <- as.vector(as.matrix(ctnames_df))
    ctnames_v <- ctnames_v[!is.na(ctnames_v)]

    deconv_results <- deconv_results[,ctnames_v]

}


rowsum(mcp_dcv[,mcp_ctnames$T_cells])

# plot correlations
