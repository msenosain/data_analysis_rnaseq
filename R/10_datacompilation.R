library(edgeR)
fc <- read.delim("counts_fcounts.txt", stringsAsFactors = FALSE)
genes <- fc[,1:6]
counts <- data.matrix(fc[,7:11])
row.names(counts) <- paste(genes$Geneid, genes$Start, genes$End, sep=".")
group <- factor( c("FL","FL","HL","HL") )
y <- DGEList(counts, genes=genes, group=group)


rnaseq_counts <- read.delim("~/Documents/Massion_lab/RNA Seq/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.count.txt")

genes <- rnaseq_counts[,1:8]
counts <- as.matrix(rnaseq_counts[,9:ncol(rnaseq_counts)])

y <- DGEList(counts, genes=genes)

# y$samples$group <- CANARY or any other classification

# y$samples$lane <- lane  (batch)

# y$samples should have all information on the samples that will be using in the analysis

###############################################################################

# First Batch
rna_3388YZ_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/RNA Seq/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.count.txt"
fpkm_3388YZ_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/RNA Seq/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.fpkm.txt"

rna_3388YZ <- read.delim(rna_3388YZ_path)
fpkm_3388YZ <- read.delim(fpkm_3388YZ_path)
# table(rna_3388YZ$Feature == fpkm_3388YZ$Feature) #TRUE
# table(colnames(rna_3388YZ) == colnames(fpkm_3388YZ)) #TRUE

p_3388YZ <- readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/RNA Seq/TMA 36 RNA DNA.xlsx", 
    sheet = "3388-YZ", col_types = c("text", 
        "text"))

x <- match(paste0('pt.', p_3388YZ$pt_ID), colnames(rna_3388YZ)[9:ncol(rna_3388YZ)]) + 8

rna_3388YZ <- rna_3388YZ[,c(1:8,x)]
fpkm_3388YZ <- fpkm_3388YZ[,c(1:8,x)]

table(colnames(rna_3388YZ)[9:70] == paste0('pt.', p_3388YZ$pt_ID)) # confirm
table(colnames(fpkm_3388YZ)[9:70] == paste0('pt.', p_3388YZ$pt_ID)) # confirm

colnames(rna_3388YZ)[9:70] = p_3388YZ$Vantage_ID
colnames(fpkm_3388YZ)[9:70] = p_3388YZ$Vantage_ID

p_3388YZ['Batch'] = 1

rna_3388YZ$Feature_gene_name1 <- NULL
fpkm_3388YZ$Feature_gene_name1 <- NULL

p_3388YZ[37, 2] <- '14301'
p_3388YZ[56, 2] <- '8356'


# Second Batch

rna_4163YZ_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/RNA Seq/RnaSeq_4163YZ_Tumor/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.count.txt"
fpkm_4163YZ_path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/RNA Seq/RnaSeq_4163YZ_Tumor/RnaSeq_3388YZ_Tumor/RnaSeq_3388YZ_Tumor.fpkm.txt"

rna_4163YZ <- read.delim(rna_4163YZ_path)
fpkm_4163YZ <- read.delim(fpkm_4163YZ_path)
# table(rna_4163YZ$Feature == fpkm_4163YZ$Feature) #TRUE
# table(colnames(rna_4163YZ) == colnames(fpkm_4163YZ)) #TRUE

p_4163YZ <- readxl::read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/RNA Seq/TMA 36 RNA DNA.xlsx", 
    sheet = "4163-YZ", col_types = c("text", 
        "text"))

x <- match(paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)), colnames(rna_4163YZ)[8:ncol(rna_4163YZ)]) + 7

rna_4163YZ <- rna_4163YZ[,c(1:7,x)]
fpkm_4163YZ <- fpkm_4163YZ[,c(1:7,x)]

table(colnames(rna_4163YZ)[8:36] == paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)))
table(colnames(fpkm_4163YZ)[8:36] == paste0('X',gsub('-', '.', p_4163YZ$Vantage_ID)))

colnames(rna_4163YZ)[8:36] = p_4163YZ$Vantage_ID
colnames(fpkm_4163YZ)[8:36] = p_4163YZ$Vantage_ID

p_4163YZ['Batch'] = 2

p_4163YZ <- p_4163YZ[-29,]
rna_4163YZ[,36] <- NULL
fpkm_4163YZ[,36] <- NULL


# Combining both 

x <- which(rna_3388YZ$Feature %in% rna_4163YZ$Feature)
x2 <- which(rna_4163YZ$Feature %in% rna_3388YZ$Feature)

rna_3388YZ <- rna_3388YZ[x,]
fpkm_3388YZ <- fpkm_3388YZ[x,]

rna_4163YZ <- rna_4163YZ[x2,]
fpkm_4163YZ <- fpkm_4163YZ[x2,]

rna_3388YZ$Feature <- factor(rna_3388YZ$Feature)
fpkm_3388YZ$Feature <- factor(fpkm_3388YZ$Feature)
rna_4163YZ$Feature <- factor(rna_4163YZ$Feature)
fpkm_4163YZ$Feature <- factor(fpkm_4163YZ$Feature)

table(rna_3388YZ$Feature == rna_4163YZ$Feature) # all should be TRUE
table(fpkm_3388YZ$Feature == fpkm_4163YZ$Feature) # all should be TRUE

rna_all <- cbind(rna_3388YZ, rna_4163YZ[,8:ncol(rna_4163YZ)])
fpkm_all <- cbind(fpkm_3388YZ, fpkm_4163YZ[,8:ncol(fpkm_4163YZ)])
p_all <- rbind(p_3388YZ, p_4163YZ)

p_all$Vantage_ID = paste0('R',gsub('-', '_', p_all$Vantage_ID))

colnames(rna_all)[8:ncol(rna_all)] = p_all$Vantage_ID
colnames(fpkm_all)[8:ncol(fpkm_all)] = p_all$Vantage_ID


# Matching with  CDE
CDE_TMA36 <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')
x <- match(p_all$pt_ID, CDE_TMA36$pt_ID) # Select only pts with RNA Seq data
pData_rnaseq <- CDE_TMA36[x,]


# Normalizing by gene length (TPM)
tpm <- function(counts, lengths) {
  x <- counts / lengths
  tpm.mat <- t( t(x) * 1e6 / colSums(x) )
  tpm.mat
}


tpm_all <- tpm(rna_all[,8:ncol(rna_all)], rna_all$Feature_length)
tpm_all_cbs <- cbind('sample_ID' = rna_all$Feature_gene_name, data.frame(tpm_all))
rownames(tpm_all) <- rna_all$Feature_gene_name

fpkm_dcv <- as.matrix(fpkm_all[,8:ncol(fpkm_all)])
fpkm_dcv_cbs <- cbind('sample_ID' = fpkm_all$Feature_gene_name, data.frame(fpkm_dcv))
rownames(fpkm_dcv) <- fpkm_all$Feature_gene_name



save(rna_all, fpkm_all, fpkm_dcv, fpkm_dcv_cbs, tpm_all, tpm_all_cbs, p_all, pData_rnaseq, file = 'rnaseq.RData')








