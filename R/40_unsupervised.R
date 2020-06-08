library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(DESeq2)

# Load data
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")

# Normalization log2(fpkm+1)
log2fpkm <- cbind(fpkm_all[,1:7], log(fpkm_all[,8:ncol(fpkm_all)] + 1, base = 2))

# Remove genes without variance between patients
variances <- apply(log2fpkm[, 8:ncol(log2fpkm)], 1, var)
sd <- apply(log2fpkm[, 8:ncol(log2fpkm)], 1, sd)
q1 <- quantile(variances, na.rm = T)["25%"]
log2fpkm <- log2fpkm[-which(is.na(variances) | variances <= q1 | sd == 0), ]


# apply combat for batch effect

# parametric adjustment
combat_edata1 = ComBat(dat=as.matrix(log2fpkm[,8:ncol(log2fpkm)]), batch=p_all$Batch, mod=NULL, par.prior=TRUE, prior.plots=TRUE)
# Batches by PCA
pca_rna <- prcomp(t(log2fpkm[,8:ncol(log2fpkm)]))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 
pca_rna <- prcomp(t(combat_edata1))
plot(pca_rna$x, col=p_all$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 


# comparison between replicates, correlation R and pval

# 11817
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_44", y = 'R4163_YZ_6',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_44", y = 'R4163_YZ_6',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 11840
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_45", y = 'R4163_YZ_28',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_45", y = 'R4163_YZ_28',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 12889
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_55", y = 'R4163_YZ_10',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_55", y = 'R4163_YZ_10',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 12890
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_6", y = 'R4163_YZ_27',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_6", y = 'R4163_YZ_27',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 12929
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_20", y = 'R4163_YZ_11',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_20", y = 'R4163_YZ_11',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 13034
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_8", y = 'R4163_YZ_12',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_8", y = 'R4163_YZ_12',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 13155
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_9", y = 'R4163_YZ_14',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_9", y = 'R4163_YZ_14',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 

# 15002
ggpubr::ggscatter(log2fpkm, x = "R3388_YZ_40", y = 'R4163_YZ_24',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 
ggpubr::ggscatter(data.frame(combat_edata1), x = "R3388_YZ_40", y = 'R4163_YZ_24',
          add = "reg.line", conf.int = TRUE, combine = TRUE,
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = 'grey')) 



# Get the average of duplicated genes
result <- WGCNA::collapseRows(combat_edata1,
                          rowGroup=log2fpkm$Feature_gene_name,
                          rowID=rownames(combat_edata1),
                          method="function",
                          methodFunction=colMeans)

data <- data.frame(result$datETcollapsed) #use this for undupervised analysis
variances <- apply(data, 1, var)
sdv <- apply(data, 1, sd)


top_unsup <- data.frame(data) %>%
  mutate(gene =rownames(.),
    sdv = sdv,
    variances = variances) %>%
  arrange(desc(variances)) %>%
  select(gene) %>%
  head(100) %>%
  pull()

filtered_data <- data.frame(data) %>%
  mutate(gene =rownames(.)) %>%
  filter(gene %in% top_unsup) %>%
  select(.,gene, everything())

rownames(filtered_data) <- filtered_data$gene

heatmap(as.matrix(filtered_data[,2:ncol(filtered_data)]), scale='row')


cormat <- Hmisc::rcorr(t(as.matrix(filtered_data[,2:ncol(filtered_data)])))

heatmap(cormat$r)

# try removing unsignificant p vals from matrix

###############################################################################

#Unsupervised/ clustering, genes and patients corplots


library(readr)
library(dplyr)
library(tidyr)

# Load RNA seq data 
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")

# Load CyTOF data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/prcnts_by_pt.RData")

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

pt_cytof <- pData_cytof$pt_ID[which(pData_cytof$pt_ID %in% unique_pt)]


p_all_cytof <- p_all %>%
  filter(., Vantage_ID %in% unique_vntg,
    pt_ID %in% pt_cytof)

fpkm_dcv_u <- fpkm_dcv[,p_all_cytof$Vantage_ID]


# Normalization log2(fpkm+1)
log2fpkm <- cbind(fpkm_all[,1:7], log( fpkm_dcv_u + 1, base = 2))

# Remove genes without variance between patients
variances <- apply(log2fpkm[, 8:ncol(log2fpkm)], 1, var)
sd <- apply(log2fpkm[, 8:ncol(log2fpkm)], 1, sd)
q1 <- quantile(variances, na.rm = T)["25%"]
log2fpkm <- log2fpkm[-which(is.na(variances) | variances <= q1 | sd == 0), ]


# apply combat for batch effect

# parametric adjustment
combat_edata1 = sva::ComBat(dat=as.matrix(log2fpkm[,8:ncol(log2fpkm)]), batch=p_all_cytof$Batch, mod=NULL, par.prior=TRUE, prior.plots=F)
# Batches by PCA
pca_rna <- prcomp(t(log2fpkm[,8:ncol(log2fpkm)]))
plot(pca_rna$x, col=p_all_cytof$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 
pca_rna <- prcomp(t(combat_edata1))
plot(pca_rna$x, col=p_all_cytof$Batch, pch=19)
legend(150, -50 , legend = c("Batch 1", "Batch 2"), pch = 19, col=c("black", "red")) 



# Get the average of duplicated genes
result <- WGCNA::collapseRows(combat_edata1,
                          rowGroup=log2fpkm$Feature_gene_name,
                          rowID=rownames(combat_edata1),
                          method="function",
                          methodFunction=colMeans)

data <- data.frame(result$datETcollapsed) #use this for undupervised analysis
variances <- apply(data, 1, var)
sdv <- apply(data, 1, sd)

# top genes with highest variances
top_unsup <- data.frame(data) %>%
  mutate(gene =rownames(.),
    sdv = sdv,
    variances = variances) %>%
  arrange(desc(variances)) %>%
  select(gene) %>%
  head(100) %>%
  pull()

filtered_data <- data.frame(data) %>%
  mutate(gene =rownames(.)) %>%
  filter(gene %in% top_unsup) %>%
  select(.,gene, everything())

rownames(filtered_data) <- filtered_data$gene
filtered_data$gene <- NULL

# Clinical vars
pData_cytof_rna <- pData_cytof
rownames(pData_cytof_rna) <- as.character(pData_cytof$pt_ID)
pData_cytof_rna <- pData_cytof_rna[p_all_cytof$pt_ID,]

ha = rowAnnotation(

    Hist_pred = as.factor(pData_cytof_rna$Hist_predominant),
    Hist_sec = as.factor(pData_cytof_rna$Hist_other_patterns),

    simple_anno_size = unit(0.5, "cm")
)

set.seed(45)
Heatmap(t(as.matrix(filtered_data)), name = "mat", row_km = 2, column_km = 2,
  heatmap_legend_param = list(color_bar = "continuous"), 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8), right_annotation = ha)



    canary = as.factor(pData_cytof_rna$CANARY), #
    age = pData_cytof_rna$Age_at_collection,
    gender = as.factor(pData_cytof_rna$Gender),
    race = as.factor(pData_cytof_rna$Race),

    ethnicity = as.factor(pData_cytof_rna$Ethnicity),
    BMI = pData_cytof_rna$BMI,
    smoking = as.factor(pData_cytof_rna$Smoking_Status),
    Pack_Years = pData_cytof_rna$Pack_Years,
    
    Age_Started = pData_cytof_rna$Age_Started,
    Age_Quit = pData_cytof_rna$Age_Quit,
    Exposure = as.factor(pData_cytof_rna$Exposure),
    prior_cancer = as.factor(pData_cytof_rna$Prior_Cancer),
    
    prior_cancer_type = as.factor(pData_cytof_rna$Prior_Cancer_Type), #*
    family_cancer = as.factor(pData_cytof_rna$Family_History_Cancer_Type),
    Chest_CT_Size = pData_cytof_rna$Chest_CT_Size,
    Chest_CT_Location = as.factor(pData_cytof_rna$Chest_CT_Location),
    
    Chest_CT_Nodule_Density = as.factor(pData_cytof_rna$Chest_CT_Nodule_Density), #
    CT_Nodule_Margination = as.factor(pData_cytof_rna$CT_Nodule_Margination),
    PET_Lesion = as.factor(pData_cytof_rna$PET_Lesion),
    fev1 = as.numeric(pData_cytof_rna$`FEV1%pred`),
    
    eighth_ed_stage = as.factor(pData_cytof_rna$`8th_ed_path_stage`),
    simplified_stage = as.factor(pData_cytof_rna$Stages_simplified),
    Path_Nodule_Size_cm = pData_cytof_rna$Path_Nodule_Size_cm,
    life_status = as.factor(pData_cytof_rna$Living_Status),
    
    Death_st = as.factor(pData_cytof_rna$Death_st),
    Recurrence_st = as.factor(pData_cytof_rna$Recurrence_st),
    Progression_st = as.factor(pData_cytof_rna$Progression_st),
    DRP_st = as.factor(pData_cytof_rna$DRP_st),

    Hist_pred = as.factor(pData_cytof_rna$Hist_predominant),
    Hist_sec = as.factor(pData_cytof_rna$Hist_other_patterns),


corrmat <- rcorr(t(as.matrix(filtered_data))) # transpose for genes


set.seed(45)
Heatmap(corrmat$r, name = "mat", row_km = 3, column_km = 3,
  heatmap_legend_param = list(color_bar = "continuous"), 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8))#, right_annotation = ha)



corrmat <- rcorr(as.matrix(filtered_data)) # transpose for genes

ha = rowAnnotation(

    Death_st = as.factor(pData_cytof_rna$Death_st),
    Recurrence_st = as.factor(pData_cytof_rna$Recurrence_st),
    Progression_st = as.factor(pData_cytof_rna$Progression_st),
    DRP_st = as.factor(pData_cytof_rna$DRP_st),

    simple_anno_size = unit(0.5, "cm")
)

set.seed(45)
Heatmap(corrmat$r, name = "mat", row_km = 2, column_km = 2,
  heatmap_legend_param = list(color_bar = "continuous"), 
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8), right_annotation = ha)



    canary = as.factor(pData_cytof_rna$CANARY), #
    age = pData_cytof_rna$Age_at_collection,
    gender = as.factor(pData_cytof_rna$Gender),
    race = as.factor(pData_cytof_rna$Race),

    ethnicity = as.factor(pData_cytof_rna$Ethnicity),
    BMI = pData_cytof_rna$BMI,
    smoking = as.factor(pData_cytof_rna$Smoking_Status),
    Pack_Years = pData_cytof_rna$Pack_Years,
    
    Age_Started = pData_cytof_rna$Age_Started,
    Age_Quit = pData_cytof_rna$Age_Quit,
    Exposure = as.factor(pData_cytof_rna$Exposure),
    prior_cancer = as.factor(pData_cytof_rna$Prior_Cancer),
    
    prior_cancer_type = as.factor(pData_cytof_rna$Prior_Cancer_Type), #*
    family_cancer = as.factor(pData_cytof_rna$Family_History_Cancer_Type),
    Chest_CT_Size = pData_cytof_rna$Chest_CT_Size,
    Chest_CT_Location = as.factor(pData_cytof_rna$Chest_CT_Location),
    
    Chest_CT_Nodule_Density = as.factor(pData_cytof_rna$Chest_CT_Nodule_Density), #
    CT_Nodule_Margination = as.factor(pData_cytof_rna$CT_Nodule_Margination),
    PET_Lesion = as.factor(pData_cytof_rna$PET_Lesion),
    fev1 = as.numeric(pData_cytof_rna$`FEV1%pred`),
    
    eighth_ed_stage = as.factor(pData_cytof_rna$`8th_ed_path_stage`),
    simplified_stage = as.factor(pData_cytof_rna$Stages_simplified),
    Path_Nodule_Size_cm = pData_cytof_rna$Path_Nodule_Size_cm,
    life_status = as.factor(pData_cytof_rna$Living_Status),
    
    Death_st = as.factor(pData_cytof_rna$Death_st),
    Recurrence_st = as.factor(pData_cytof_rna$Recurrence_st),
    Progression_st = as.factor(pData_cytof_rna$Progression_st),
    DRP_st = as.factor(pData_cytof_rna$DRP_st),

    Hist_pred = as.factor(pData_cytof_rna$Hist_predominant),
    Hist_sec = as.factor(pData_cytof_rna$Hist_other_patterns),





