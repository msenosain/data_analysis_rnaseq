environment_set <- function(){
    library(EnhancedVolcano)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(DESeq2)
    library(sva)
    library(limma)
    library(ComplexHeatmap)
    library(fgsea)
    library(tibble)
    library(gage)
    library(ggplot2)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    # Check BiocManager::valid()
}

preprocess_rna <- function(path_rnaseq, 
    correct_batch=TRUE, 
    lowvargenesrm = TRUE,
    xychr_rm = TRUE,
    prot_coding_only = FALSE){
    # Load data
    load(path_rnaseq)

    # Remove low quality samples
    qc_R4163 <- paste0('R4163_YZ_', c('4','12','13','14','15','27'))
    qc_R3388 <- paste0('R3388_YZ_', c('8', '45'))

    # Remove duplicates
    dupl <- c('11817', '12889', '12929', '15002') #11840 13034 12890 13155
    unique_vntg <- p_all %>%
      filter(., !(pt_ID %in% dupl &
        grepl("R4163",Vantage_ID)),
        !(Vantage_ID %in% c(qc_R4163,qc_R3388))) %>%
      dplyr::select(., Vantage_ID) %>%
      pull()

    p_all <- p_all[which(p_all$Vantage_ID %in% unique_vntg),]
    rna_all <- cbind(rna_all[,1:7], rna_all[,unique_vntg])
    pData_rnaseq <- pData_rnaseq[!duplicated(pData_rnaseq$pt_ID),]
    pData_rnaseq$pt_ID <- as.character(pData_rnaseq$pt_ID)
    pData_rnaseq <- pData_rnaseq[match(p_all$pt_ID, pData_rnaseq$pt_ID),]

    # Remove low variance genes from counts
    if(lowvargenesrm){
        #variances <- apply(rna_all[, 8:ncol(rna_all)], 1, var)
        #sdv <- apply(rna_all[, 8:ncol(rna_all)], 1, sd)
        #q1 <- quantile(variances, na.rm = T)["25%"]
        #idx <- which(is.na(variances) | variances <= q1 | sdv == 0)
        x<- rna_all[,8:ncol(rna_all)]
        idx <- edgeR::filterByExpr(x)
        rna_all <- rna_all[idx, ]
    }

    if(xychr_rm){
        idx <- which(rna_all$Feature_chr %in% c('chrX', 'chrY'))
        rna_all <- rna_all[-idx, ]
    }

    if(prot_coding_only){
        idx <- which(rna_all$Feature_gene_biotype == 'protein_coding')
        rna_all <- rna_all[idx, ]
    }    
    
    counts_all <- rna_all[, 8:ncol(rna_all)]
    rownames(counts_all) <- rna_all$Feature

    # Normalize
    p_all$Batch <- as.factor(p_all$Batch)
    dds <- DESeqDataSetFromMatrix(
      countData = counts_all,
      colData = p_all,
      design = ~Batch, tidy = F  
    )
    dds <- estimateSizeFactors(dds)
    sizeFactors(dds)
    vsd <- vst(dds)
    vsd_mat <- assay(vsd)

    # Assessing batch effect
    pbatch_bf <- plotPCA(vsd, "Batch") + labs(fill = "Batch") + ggtitle("Batch raw")

    # Save pre-processed data in a list
    dge_preprocessed <- list(p_all=p_all, rna_all=rna_all, 
        pData_rnaseq=pData_rnaseq, counts_all=counts_all, vsd_mat=vsd_mat,
        pbatch_bf=pbatch_bf)

    # Remove batch effect
    if(correct_batch){
        assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$Batch)
        pbatch_af <- plotPCA(vsd, "Batch") + labs(fill = "Batch") + ggtitle("Batch after BE removal")
        vsd_mat <- assay(vsd)
        dge_preprocessed[['vsd_mat']] <- vsd_mat
        dge_preprocessed[['pbatch_af']] <- pbatch_af
    }
    
    dge_preprocessed
}

label_by_gene <- function(normalized_data, gene, metadata, extremes_only=FALSE){
    if(extremes_only){
        data_t <- data.frame(t(normalized_data))
        first_q <- quantile(data_t[,gene])[["25%"]]
        third_q <- quantile(data_t[,gene])[["75%"]]
        data_t$Condition <- "normal"
        data_t$Condition[data_t[,gene] < first_q] <- "low"
        data_t$Condition[data_t[,gene] > third_q] <- "high"

    } else {
        data_t <- data.frame(t(normalized_data))
        mid_q <- quantile(data_t[,gene])[["50%"]]
        data_t$Condition <- "low"
        data_t$Condition[data_t[,gene]  > mid_q] <- "high"

    }

    meta_data <- cbind(metadata, Condition=data_t$Condition)
    meta_data$Batch <- factor(meta_data$Batch)
    meta_data$Gender <- factor(meta_data$Gender)
    meta_data$Condition <- factor(meta_data$Condition)

    meta_data
}

DE_analysis <- function(ls_preprocessed, 
    GeneBased=TRUE, 
    pDataBased=FALSE,
    NewCondition=FALSE,
    NewCondition_df, # condition label into p_all
    cond_nm, # gene or column name
    two_levels=c('high','low'), # high low, dead alive...
    reference, # low, alive
    extremes_only=TRUE
    ){

    # Unlist
    p_all=ls_preprocessed$p_all
    rna_all=ls_preprocessed$rna_all
    pData_rnaseq=ls_preprocessed$pData_rnaseq
    counts_all=ls_preprocessed$counts_all
    vsd_mat=ls_preprocessed$vsd_mat

    message('Unlist done')


    if(GeneBased){
        meta_data <- label_by_gene(vsd_mat, gene=cond_nm, p_all, extremes_only=extremes_only)
        meta_data <- meta_data %>% filter(Condition != "normal")
    }
    
    if(pDataBased){
        meta_data <- pData_rnaseq %>% dplyr::select(pt_ID, all_of(cond_nm))
        meta_data$pt_ID <- as.character(meta_data$pt_ID)
        p_all$pt_ID <- as.character(p_all$pt_ID)
        meta_data <- inner_join(p_all, meta_data, by='pt_ID')
        names(meta_data)[names(meta_data) == cond_nm] <- 'Condition'
        meta_data <- meta_data %>% filter(Condition %in% two_levels)
    }

    if(NewCondition){
        meta_data <- NewCondition_df
        print(meta_data)
        names(meta_data)[names(meta_data) == cond_nm] <- 'Condition'
        print(meta_data)
        meta_data$Condition <- as.character(meta_data$Condition)
        print(meta_data)
        meta_data <- meta_data %>% filter(Condition %in% two_levels)
        print(meta_data)
    }

    message('Labeling done')

    counts_all <- as.matrix(counts_all[,meta_data$Vantage_ID])
    pData_rnaseq <- pData_rnaseq %>% filter(pt_ID %in% meta_data$pt_ID)
    vsd_mat <- vsd_mat[,meta_data$Vantage_ID]

    message('Filtering done')

    dds <- DESeqDataSetFromMatrix(
      countData = counts_all,
      colData = meta_data,
      design = ~Condition, tidy = F) #design = ~Batch + Condition, tidy = F)

    message('Design done')

    dds <- estimateSizeFactors(dds)

    # Generate a transformed matrix with gene symbols
    ens2symbol <- data.frame(cbind(ENSEMBL=as.character(rna_all$Feature), 
        symbol=as.character(rna_all$Feature_gene_name)))
    vsd_mat_sym <- data.frame(vsd_mat) %>% mutate(gene=rownames(.)) %>% as_tibble()
    vsd_mat_sym <- inner_join(vsd_mat_sym, ens2symbol, by=c("gene"="ENSEMBL"))
    rownames(vsd_mat_sym) <- make.names(vsd_mat_sym$symbol, unique=TRUE)

    message('vsd symbols done')

    dds$Condition <- relevel(dds$Condition, ref = reference)
    dds <- DESeq(dds, parallel = F)

    message('DESeq done')

    res <- results(dds, contrast = c('Condition', two_levels))
    #res <- lfcShrink(dds, contrast = c('Condition',two_levels), res=res, type = 'normal')
    res_df <- data.frame(res) %>% mutate(gene=rownames(.)) %>% as_tibble()

    res_df <- res_df %>% 
      inner_join(., ens2symbol, by=c("gene"="ENSEMBL")) %>%
      data.frame()
    rownames(res_df) <- make.names(res_df$symbol, unique=TRUE)

    message('res symbols done')

    DE_res <- list(dds=dds, vsd_mat=vsd_mat, vsd_mat_sym=vsd_mat_sym, 
        ens2symbol=ens2symbol, res=res, res_df=res_df, 
        pData_rnaseq=pData_rnaseq, meta_data=meta_data)

    message('list done')

    DE_res
}

fgsea_analysis <- function(DE_res){

    # Unlist
    res_df=DE_res$res_df

    ranks <- res_df %>%
      dplyr::select(symbol,stat) %>%
      na.omit() %>% 
      distinct() %>% 
      group_by(symbol) %>% 
      summarize(stat=mean(stat)) %>%
      arrange(stat)

    ranks <- deframe(ranks)

    hallmark_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/h.all.v7.1.symbols.gmt'
    c1_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c1.all.v7.1.symbols.gmt'
    c2_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c2.all.v7.1.symbols.gmt'
    c3_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c3.all.v7.1.symbols.gmt'
    c4_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c4.all.v7.1.symbols.gmt'
    c5_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c5.all.v7.1.symbols.gmt'
    c6_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c6.all.v7.1.symbols.gmt'
    c7_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/c7.all.v7.1.symbols.gmt'
    msig_path <- '/Users/senosam/Documents/Massion_lab/RNASeq_summary/GSEA/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt'

    # For REACTOME pathways (only work with ENTREZID)

    fgsea_fixed <- function(pthw_path, reactome=FALSE, ranks=ranks){
        if(reactome){
            library(org.Hs.eg.db)
            hs <- org.Hs.eg.db
            my.symbols <- c(names(ranks))
            entrz <- AnnotationDbi::select(hs, 
                keys = my.symbols,
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL")
            entrz <- na.omit(entrz)
            ranks <- ranks[entrz$SYMBOL]
            names(ranks) <- entrz$ENTREZID
            pth <- reactomePathways(names(ranks))
        }else{
            pth <- gmtPathways(pthw_path)
        }
        res <- fgsea(pathways=pth, stats=ranks, eps = 0) %>% #multilevel, add nPermSimple = 10000
                    as_tibble %>%
                    arrange(padj, desc(abs(NES)))
        res$state <- ifelse(res$NES > 0, "up", "down")
        res$leadingEdge <- sapply(res$leadingEdge, . %>% {
          str_c(., collapse = " ")})

        res
    }



    res_hm <- fgsea_fixed(pthw_path=hallmark_path, reactome=FALSE, ranks=ranks)
    res_c1 <- fgsea_fixed(pthw_path=c1_path, reactome=FALSE, ranks=ranks)
    res_c2 <- fgsea_fixed(pthw_path=c2_path, reactome=FALSE, ranks=ranks)
    res_c3 <- fgsea_fixed(pthw_path=c3_path, reactome=FALSE, ranks=ranks)
    res_c4 <- fgsea_fixed(pthw_path=c4_path, reactome=FALSE, ranks=ranks)
    res_c5 <- fgsea_fixed(pthw_path=c5_path, reactome=FALSE, ranks=ranks)
    res_c6 <- fgsea_fixed(pthw_path=c6_path, reactome=FALSE, ranks=ranks)
    res_c7 <- fgsea_fixed(pthw_path=c7_path, reactome=FALSE, ranks=ranks)
    res_msg <- fgsea_fixed(pthw_path=msig_path, reactome=FALSE, ranks=ranks)
    res_rtm <- fgsea_fixed(reactome=TRUE, ranks=ranks)

    fgsea_res <- list(res_hm=res_hm, res_c1=res_c1, res_c2=res_c2, 
        res_c3=res_c3, res_c4=res_c4, res_c5=res_c5, res_c6=res_c6,
        res_c7=res_c7, res_msg=res_msg, res_rtm=res_rtm)

    fgsea_res
}

##########################################################################
# Pathway analysis with KEGG pathways and GO
##########################################################################
kegg_go <- function(DE_res, kegg = TRUE, GO = FALSE){
    
    library(gage)
    library(pathview)
    library(gageData)
    data(kegg.sets.hs)
    data(go.sets.hs)

    # Function to put results in a df
    ls2df <- function(pth_res){
        pth_res <- pth_res[['greater']] %>%
                as.data.frame() %>%
                na.omit()
        pth_res$state <- ifelse(pth_res$stat.mean > 0, "up", "down")
        pth_res
    }  

    # setting up the input
    res_df=DE_res$res_df
    res_df$ENSEMBL <- sapply(strsplit(res_df$gene, "\\."), "[[", 1)

    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    ensmbl <- res_df$ENSEMBL
    entrz <- AnnotationDbi::select(hs, 
        keys = ensmbl,
        columns = c("ENTREZID", "ENSEMBL"),
        keytype = "ENSEMBL") %>% na.omit()
    entrz <- merge(res_df, entrz, by='ENSEMBL')

    fc <- entrz$log2FoldChange
    names(fc) <- entrz$ENTREZID

    # running the KEGG pathway analysis
    if(kegg){
        keggres <- gage(fc, gsets=kegg.sets.hs, same.dir=FALSE)
        keggres <- ls2df(keggres)
        return(keggres)
    }

    # Running GO 
    if(GO){
        gores = gage(fc, gsets=go.sets.hs, same.dir=FALSE)
        gores <- ls2df(gores) 
        return(gores)       
    }
    #https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
}

# plots
heatmap_200 <- function(res_df, vsd_mat, meta_data, pData_rnaseq, n_genes=200,
    pval_cutoff = 0.05, l2fc_cutoff=1.5, ha_custom=NULL, row_km=2, scale_mat = F){
    res_df <- data.frame(res_df) %>%
        na.omit() %>%
        filter(abs(log2FoldChange) > l2fc_cutoff) %>%
        filter(pvalue < pval_cutoff) %>%
        arrange(pvalue)

    if(nrow(res_df)>n_genes){
        gene_list <- res_df %>%
            head(n_genes) %>%
            dplyr::select(gene) %>%
            pull()
    }else{
        gene_list <- res_df %>%
            dplyr::select(gene) %>%
            pull()        
    }

    # Use normalized data (vsd_mat) for plot
    filtered_res <- data.frame(vsd_mat) %>%
      filter(gene %in% gene_list) %>%
      as.data.frame()
    rownames(filtered_res) <- make.names(filtered_res$symbol, unique=TRUE)
    filtered_res$gene <- NULL
    filtered_res$symbol <- NULL

    ha = HeatmapAnnotation(

        CANARY = as.factor(pData_rnaseq$CANARY),
        #gender = as.factor(pData_rnaseq$Gender),
        Condition = as.factor(meta_data$Condition),

        simple_anno_size = unit(0.5, "cm")
    )

    if(!is.null(ha_custom)){
        ha <- ha_custom
    }

    if(scale_mat){
        filtered_res <- t(scale(t(as.matrix(filtered_res))))
    }
    print(Heatmap(filtered_res, name = "mat", 
        #column_km = 2, 
      #row_km = row_km,
      column_split =as.factor(meta_data$Condition),
      heatmap_legend_param = list(color_bar = "continuous"), 
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8), top_annotation = ha))
}

volcano_plot <- function(res_df, gene=NULL, p_title, pCutoff=0.05, FCcutoff=1.5){

    if (!is.null(gene)){
        k <- which(res_df$gene %in% gene)
        res_df <- res_df[-k,]
    }
    res_df <- na.omit(res_df) #remove NA, these are outliers
    print(EnhancedVolcano(res_df,
        lab = rownames(res_df),
        x = 'log2FoldChange',
        y = 'pvalue',
        pCutoff=pCutoff,
        FCcutoff=FCcutoff,
        xlim = c(-5,5),
        pointSize=2,
        labSize=4,
        title = p_title))
}

fgsea_plot <- function(fgsea_res, pathways_title, cutoff = 0.05, 
    max_pathways = 30, condition_name){

        color_levels <- function(fgsea_res) {
            colors <- c()
            if (any(fgsea_res$state == "down")) {
              colors <- c(colors, "lightblue")
            }
            if (any(fgsea_res$state == "up")) {
              colors <- c(colors, "#DC143C")
            }
            colors
        }

        # Add * code for p vals
        fgsea_res$pvlabel <- '*'
        fgsea_res$pvlabel[which(fgsea_res$padj <0.01 & fgsea_res$padj>0.001)] <- '**'
        fgsea_res$pvlabel[which(fgsea_res$padj<0.001)] <- '***'

        if (!is.null(cutoff)) {
            fgsea_res <- fgsea_res %>% filter(padj < cutoff)
        }
        
        curated_pathways <- fgsea_res %>%
                arrange(desc(abs(NES))) %>%
                dplyr::slice(1:max_pathways)
        curated_pathways['leadingEdge'] <- NULL
        print(ggplot(curated_pathways, aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill = state), width = 0.5, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            geom_label(aes(label = pvlabel), size = 3, alpha = 0.75) +
            coord_flip() +
            labs(
                x = 'Pathway', 
                y = "Normalized Enrichment Score",
                title = str_c(pathways_title, " pathways: ", condition_name),
                subtitle = str_c("(Cutoff: p.adj <", cutoff, ")")
            ) +
            theme_bw() +
            scale_fill_manual(values = color_levels(curated_pathways)))

        fgsea_res <- fgsea_res %>% 
                #dplyr::select(-leadingEdge, -ES) %>% 
                arrange(desc(abs(NES)))
                #DT::datatable())

        knitr::kable(fgsea_res)
}

keggGO_plot <- function(keggGO_res, pathways_title, cutoff = 0.05, 
    max_pathways = 30, condition_name, pval_colnm){

        color_levels <- function(fgsea_res) {
            colors <- c()
            if (any(keggGO_res$state == "down")) {
              colors <- c(colors, "lightblue")
            }
            if (any(keggGO_res$state == "up")) {
              colors <- c(colors, "#DC143C")
            }
            colors
        }

        # Add * code for p vals
        keggGO_res$pvlabel <- '*'
        keggGO_res$pvlabel[which(keggGO_res[,pval_colnm] <0.01 & keggGO_res[,pval_colnm]>0.001)] <- '**'
        keggGO_res$pvlabel[which(keggGO_res[,pval_colnm]<0.001)] <- '***'

        # Add column for pathway and pathway_id
        pathway_id <- sapply(strsplit(rownames(keggGO_res), " "), "[[", 1)
        pathway<- sub("^(?:\\S+\\s+)", "\\1", rownames(keggGO_res), perl = TRUE)
        keggGO_res <- cbind(pathway_id, pathway, keggGO_res)
        rownames(keggGO_res)<-NULL

        if (!is.null(cutoff)) {
            keggGO_res$pvselect <- keggGO_res[,pval_colnm]
            keggGO_res <- keggGO_res %>% filter(pvselect < cutoff)
        }
        
        curated_pathways <- keggGO_res %>%
                arrange(desc(abs(stat.mean))) %>%
                dplyr::slice(1:max_pathways)

        print(ggplot(curated_pathways, aes(reorder(pathway, stat.mean), stat.mean)) +
            geom_col(aes(fill = state), width = 0.5, color = "black") +
            scale_size_manual(values = c(0, 1), guide = "none") +
            geom_label(aes(label = pvlabel), size = 3, alpha = 0.75) +
            coord_flip() +
            labs(
                x = '', 
                y = "stats.mean",
                title = str_c(pathways_title, ": ", condition_name),
                subtitle = str_c("(Cutoff: ", pval_colnm, " < ", cutoff, ")")
            ) +
            theme_bw() +
            scale_fill_manual(values = color_levels(curated_pathways)))

        keggGO_res <- keggGO_res %>% 
                dplyr::select(-exp1, -pvselect) %>% 
                arrange(desc(abs(stat.mean)))
                #DT::datatable())

        knitr::kable(keggGO_res)
}


#x <- kegg_go(DE_res, kegg = TRUE, GO = FALSE)
#keggGO_plot(x, pathways_title='KEGG', cutoff = 0.05, 
#                        max_pathways = 30, condition_name='I VS A', pval_colnm='p.val')


