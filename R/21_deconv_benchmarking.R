library(readr)
library(dplyr)
library(tidyr)
source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')

# Load deconvolution data rna_cytof

# TPM
mcp_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/MCP/mcp_dcv.txt", row.names=1)))
qts_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/QTS/qts_dcv.txt", row.names=1)))
cbs_dcv <- read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/CBS/CIBERSORT.rna_cytof_tpm.txt", row.names=1)
xcell_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/XCELL/xCell_rnaseq_tpm_xCell_1211060320.txt", row.names=1)))

# FPKM
mcp_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/MCP/mcp_fpkm_dcv.txt", row.names=1)))
qts_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/QTS/qts_fpkm_dcv.txt", row.names=1)))
cbs_dcv <- read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/CBS/CIBERSORT.rna_cytof_fpkm.txt", row.names=1)
xcell_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_cytof/XCELL/xCell_rnaseq_fpkm_xCell_1151060320.txt", row.names=1)))


# Load RNA seq data 
ls_preprocessed <- preprocess_rna(path_rnaseq = 'rnaseq.RData', correct_batch = T, correct_gender = T)

# Load CyTOF data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/prcnts_by_pt.RData")

prcnt_pt <- prcnt_pt %>%
  tibble::rownames_to_column('pt_ID') %>%
  mutate(., Epithelial = rowSums(prcnt_pt[,9:ncol(prcnt_pt)])) %>%
  dplyr::select(., Epithelial, everything()) %>%
  dplyr::select(., -colnames(prcnt_pt)[8:ncol(prcnt_pt)]) %>%
  tibble::column_to_rownames('pt_ID')

unique_pt <- ls_preprocessed$p_all %>%
  filter(., Vantage_ID %in% rownames(mcp_dcv)) %>%
  dplyr::select(., pt_ID) %>%
  pull()

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
prcnt_pt['Non_immune'] <- rowSums(prcnt_pt[,1:3])
prcnt_pt['method'] <- 'cytof'


# Load excel cell types names
mcp_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "mcp")
qts_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "qts")
cbs_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "cbs")
xcell_ctnames <- readxl::read_excel("/Users/senosam/Documents/Massion_lab/RNASeq_summary/deconvolution/deconv_ctnames.xlsx", 
    sheet = "xcell")

mcp_ctnames[]<- lapply(mcp_ctnames, function(x) gsub('\\+|\\-|\\s|\\(|\\)', '.', x))
qts_ctnames[]<- lapply(qts_ctnames, function(x) gsub('\\+|\\-|\\s|\\(|\\)', '.', x))
cbs_ctnames[]<- lapply(cbs_ctnames, function(x) gsub('\\+|\\-|\\s|\\(|\\)', '.', x))
xcell_ctnames[]<- lapply(xcell_ctnames, function(x) gsub('\\+|\\-|\\s|\\(|\\)', '.', x))

# match names
rearrange_dcv <- function(ctnames_df, deconv_results, method_name){
    ctnames_v <- as.vector(as.matrix(ctnames_df))
    ctnames_v <- ctnames_v[!is.na(ctnames_v)]

    deconv_results <- deconv_results[,ctnames_v]

    new_dcv <- data.frame(matrix(0, 
        nrow = nrow(deconv_results), 
        ncol = ncol(ctnames_df), 
        dimnames = list(rownames(deconv_results), colnames(ctnames_df))))

    for(i in colnames(new_dcv)){
        v <- as.vector(as.matrix(ctnames_df[,i]))
        v <- v[!is.na(v)]
        if (length(v)==1){
            new_dcv[,i] <- deconv_results[,v]
        } else if (length(v)>1) {
            new_dcv[,i] <- rowSums(deconv_results[,v])
        } else if(length(v)==0) {
            new_dcv[,i] = NA
        }
    }
    new_dcv['method'] <- method_name
    new_dcv
}

mcp_new <- rearrange_dcv(mcp_ctnames, mcp_dcv, 'MCP')
qts_new <- rearrange_dcv(qts_ctnames, qts_dcv, 'QTS')
cbs_new <- rearrange_dcv(cbs_ctnames, cbs_dcv, 'CBS')
xcell_new <- rearrange_dcv(xcell_ctnames, xcell_dcv, 'XCL')


# plot correlations


all_dt <- rbind(mcp_new, qts_new, cbs_new, xcell_new)

all_dt <- reshape::melt(all_dt)

all_dt <- all_dt[order(all_dt$method),]

ref <- reshape::melt(prcnt_pt)

ref <- rbind(ref,ref,ref,ref)

all_dt <- cbind(all_dt, 'ref'= ref$value)


ggpubr::ggscatter(all_dt, x = "ref", y = "value",
          add = "reg.line", conf.int = F, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'blue'),
          xlab = "CyTOF fraction", ylab = 'Estimated fraction') +
    ggplot2::facet_grid(vars(method), vars(variable), scales = "free")

