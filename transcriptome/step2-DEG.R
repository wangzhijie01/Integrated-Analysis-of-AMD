rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load("step1_count_anno.Rdata")


###AMD vs CONTROL macular RPE
##黄斑区
source('./functions.R')
if(F){
  exprSet.count=counts
  pd=pdata
  pd=pd[grep("RPE, Macula",pd$tissue),]
  dim(pd)
  table(pd$tissue,pd$group)
  group_list=pd$group
  exprSet=exprSet.count[,pd$geo_accession]
  dim(exprSet)
  suppressPackageStartupMessages(library(edgeR))
  d <- DGEList(counts=exprSet,group=factor(group_list))
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  dge=d
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(1,-1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  write.csv(nrDEG,"setp2_nrdeg_edger.csv")
  head(nrDEG)
  edgeR_DEG =nrDEG 
  nrDEG=edgeR_DEG[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'edgeR-macular',group_list=group_list,logFC_cutoff=log2(1.5))
  save(exprSet,nrDEG,pd,file = "step2-DEG-macular.Rdata")
}



##非黄斑区
if(F){
  pd=pdata
  pd=pd[grep("RPE, non-Macula",pd$tissue),]
  dim(pd)
  table(pd$tissue,pd$group)
  group_list=pd$group
  exprSet=exprSet.count[,pd$geo_accession]
  dim(exprSet)
  suppressPackageStartupMessages(library(edgeR))
  d <- DGEList(counts=exprSet,group=factor(group_list))
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  dge=d
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(1,-1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  edgeR_DEG =nrDEG 
  nrDEG=edgeR_DEG[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'edgeR-non macular',group_list=group_list,logFC_cutoff=log2(1.5))
  save(exprSet,nrDEG,pd,file = "step2-DEG-non macular.Rdata")
}
