
# Reactome pathways

reapaths <- read.table('/Users/senosam/Documents/Repositories/Research/DH_project01/data/Reactome2020/NCBI2Reactome_All_Levels.txt')
genetopaths <- read.delim('/Users/senosam/Documents/Repositories/Research/DH_project01/data/Reactome2020/Ensembl2Reactome_PE_Pathway.txt', header=F, stringsAsFactors = F)

dt <- ls_preprocessed$vsd_mat
rownames(dt) <- sapply(strsplit(rownames(dt), "\\."), "[[", 1)
DE_genes <- 

genetopaths=genetopaths[which(genetopaths$V8=='Homo sapiens'),]
genetopaths$V3=sapply(genetopaths$V3,function(x){return(unlist(strsplit(x, split=' \\['))[1])})

allgenes=rownames(dt)
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
allgenes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=allgenes, mart= ensembl)

##summarise expression per pathway
paths=unique(genetopaths$V6)

pathexpl=list()

for (i in 1:length(paths)){
  #print(i)
  ens=unique(genetopaths$V1[which(genetopaths$V6==paths[i])])
  genname=unique(allgenes.with.id[which(allgenes.with.id$ensembl_gene_id %in% ens),1])
  genname_com=intersect(genname, rownames(dt))
  if (length(genname_com)>1){
  pathexpl[[paths[i]]]=colMeans(dt[genname_com,])
  }
  if (length(genname_com)<1){
    pathexpl[[paths[i]]]=NA
  } 
  if (length(genname_com)==1){
    pathexpl[[paths[i]]]=dt[genname_com,]
  } 
  
}

pathexpdf=Reduce(rbind, pathexpl)
rownames(pathexpdf)=names(pathexpl)

pathexpdf=pathexpdf[complete.cases(pathexpdf),]


pheatmap(cor(pathexpdf), show_colnames = T, show_rownames=F)



pathexpdf <- pathexpdf %>%
  data.frame(.) %>%
  mutate(OE_ave = (OverExpression1+OverExpression2+OverExpression3)/3,
    WT_ave = (WT1 + WT2 + WT3)/3)


# REACTOME vignette

library(ReactomePA)
library(topGO)
data(geneList)
de <- names(geneList)[abs(geneList) > 0.5]
head(de)

library(reactome.db)

res <- data.frame(DE_res$res)
rownames(res) <- rownames(dt) <- sapply(strsplit(rownames(res), "\\."), "[[", 1)
res2 <- res[ rownames(res) %in% keys( reactome.db, "ENTREZID" ) & !is.na( res$padj ) , ]
head(res2)

