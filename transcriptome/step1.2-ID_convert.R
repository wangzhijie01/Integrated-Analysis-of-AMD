rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load('step1-output.Rdata')

# 初步过滤
dim(exprSet.count)
counts=na.omit(exprSet.count)
counts=counts[apply(counts, 1,function(x){sum(x>1)>10}),]
dim(counts)#[1] 33006   266
counts[1:4,1:4]

if (F) {
  #ID转换
  #去除空的ENTREZ
  library('org.Hs.eg.db')
  keytypes(org.Hs.eg.db)
  gene_map<-mapIds(org.Hs.eg.db, as.character(rownames(counts)), 'SYMBOL', keytype ='ENSEMBL',multiVals = "first")
  gene_map=data.frame(ENSEMBL=names(gene_map),SYMBOL=as.character(gene_map))
  
  dat =counts
  dat$median=apply(dat,1,median)
  dat$ENSEMBL = rownames(dat)
  dat = merge(dat,gene_map,by="ENSEMBL")
  dat = na.omit(dat)
  
  dat=dat[order(dat$SYMBOL,dat$median,decreasing = T),]#把dat$symbol按照dat$median排序
  dat=dat[!duplicated(dat$SYMBOL),]#取出不重复的dat$symbol
  rownames(dat)=dat$SYMBOL
  dat=dat[,-c(1,ncol(dat)-1,ncol(dat))]
  dim(dat) 
  counts=dat
}

save(counts,pdata,file="step1_count_anno.Rdata")
