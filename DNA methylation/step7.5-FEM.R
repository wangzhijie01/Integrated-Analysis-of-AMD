

library(FEM)
options(stringsAsFactors = F)
load(file = "step2-champ_myNorm.Rdata")##载入beta矩阵
load("D:/bioinformation/R/甲基化/AMD/转录组/step1_count_anno.Rdata")##counts矩阵
load("D:/bioinformation/R/甲基化/AMD/转录组/step1-output.Rdata")##转录组表达量归一化后的矩阵
data(hprdAsigH) ##为了获取adj矩阵


#count归一化-edegR
if (F) {
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
  d <- calcNormFactors(d, method = "TMM")
  tmm <- cpm(d,log = T)
  library(org.Hs.eg.db)
  gene_map<-mapIds(org.Hs.eg.db, as.character(rownames(tmm)), 'ENTREZID', keytype ='ENSEMBL',multiVals = "first")
  rownames(tmm)=gene_map
  tmm=tmm[!is.na(rownames(tmm)),]
  # 筛选macular部分
  {
    pd=pdata[grep("RPE, Macula",pdata$tissue),]
    dim(pd)
    table(pd$tissue,pd$group)
    tmm=tmm[,pd$geo_accession]
    dim(tmm)
    exp.Normaliza=tmm
  }
}

# 直接提取RPKM值
if (F) {
  exprSet  =exprSet.RPKM
  pd=pdata
  pd=pd[grep("RPE, Macula",pd$tissue),]
  dim(pd)
  table(pd$tissue,pd$group)
  group_list=pd$group
  exprSet=exprSet.RPKM[,pd$geo_accession]
  keep <- rowSums(exprSet>0)>=ncol(exprSet)/2
  table(keep)
  exprSet=exprSet[keep,]
  dim(exprSet)
  library(org.Hs.eg.db)
  gene_map<-mapIds(org.Hs.eg.db, as.character(rownames(exprSet)), 'ENTREZID', keytype ='ENSEMBL',multiVals = "first")
  gene_map=data.frame(ENSEMBL=names(gene_map),ENTREZID=gene_map)
  gene_map$median=apply(exprSet,1,median)
  gene_map=gene_map[order(gene_map$ENTREZID,gene_map$median,decreasing = T),]
  gene_map=gene_map[!duplicated(gene_map$ENTREZID),]
  gene_map=na.omit(gene_map)
  exprSet=exprSet[gene_map$ENSEMBL,]
  exprSet[1:4,1:4]
  rownames(exprSet)=gene_map$ENTREZID
  exp.Normaliza=exprSet
}

##甲基化
dnaM.m=myNorm2
phenoM.v=myLoad$pd$Sample_Group
adj.m=hprdAsigH.m
statM.o <- GenStatM(dnaM.m,phenoM.v,"450k");
statM.o[["cont"]]
##要转变方向
statM.o[["top"]][[1]]$logFC=statM.o[["top"]][[1]]$logFC*-1
statM.o[["top"]][[1]]$t=statM.o[["top"]][[1]]$t*-1

##转录组
exp.m=exp.Normaliza
pheno.v=pd$group
statR.o <- GenStatR(exp.m,pheno.v)
statR.o[["cont"]]
##要转变方向
statR.o[["top"]][[1]]$logFC=statR.o[["top"]][[1]]$logFC*-1
statR.o[["top"]][[1]]$t=statR.o[["top"]][[1]]$t*-1

# We then need to check which contrast/column in statM.o$cont and statR.o$cont corresponds to the one of
# interest, and assign this integer to cM and cR, respectivel
DoIntFEM450k.o=DoIntFEM450k(statM.o,statR.o,adj.m,cM=1,cR=1)
# intFEM.o <- list(statM=DoIntFEM450k.o$statM,statR=DoIntFEM450k.o$statR,adj=DoIntFEM450k.o$adj)

DoFEMtoy.o <- DoFEMbi(DoIntFEM450k.o,nseeds=100,gamma=0.5,nMC=1000,
                      sizeR.v=c(1,100),minsizeOUT=10,writeOUT=TRUE,nameSTUDY="AMD-Control",ew.v=NULL);

save(DoFEMtoy.o,file="FEM_result.Rdata")
# Use the following command to display their sizeand elements.
DoFEMtoy.o$fem

# The details of the modules can also be seen using:
DoFEMtoy.o$topmod


library("marray");
library("corrplot");
setwd("./FEM/")
# 将所有图片生成导出
for(m in 1:length(names(DoFEMtoy.o$topmod))){
  a=DoFEMtoy.o$topmod[[m]]
  # a=a[!a$`stat(Int)`== 0,]
  # a=a[1:ceiling(nrow(a)/2),]#只取一半
  a=a[nrow(a)-3:nrow(a),]
  FemModShow(a,  name=names(DoFEMtoy.o$topmod)[m],DoFEMtoy.o)
  }

for(m in 1:length(names(DoFEMtoy.o$topmod))){
  a=DoFEMtoy.o$topmod[[m]]
  # a=a[!a$`stat(Int)`== 0,]
  FemModShow(a,
             name=paste0(names(DoFEMtoy.o$topmod)[m],"_FULL"),DoFEMtoy.o)}
setwd("..")
