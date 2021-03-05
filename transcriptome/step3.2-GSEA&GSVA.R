rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)


library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

load("step1_count_anno.Rdata")


# 提取macular样本
{
  exprSet.count=counts
  pd=pdata
  pd=pd[grep("RPE, Macula",pd$tissue),]
  dim(pd)
  table(pd$tissue,pd$group)
  group_list=pd$group
  exprSet=exprSet.count[,pd$geo_accession]
  dim(exprSet)
  exprSet.cpm=edgeR::cpm(exprSet)
}
#载入DEG数据
load(file = "step2-DEG-macular.Rdata")



# GSEA
{
  library(org.Hs.eg.db)
  gene_map<-mapIds(org.Hs.eg.db, as.character(rownames(nrDEG)), 'SYMBOL', keytype ='ENSEMBL',multiVals = "first")
  nrDEG$SYMBOL=gene_map
  nrDEG=na.omit(nrDEG)
  nrDEG=nrDEG[!duplicated(nrDEG$SYMBOL),]
  
  geneList=nrDEG$log2FoldChange
  names(geneList)=nrDEG$SYMBOL
  geneList=sort(geneList,decreasing = T)
  #选择gmt文件（MigDB中的全部基因集）
  d='D:/bioinformation/MsigDB/symbols/'
  gmts <- list.files(d,pattern = 'all')
  gmts
  gmts=c("h.all.v7.1.symbols.gmt","c2.cp.kegg.v7.1.symbols.gmt" )
  #GSEA分析
  library(GSEABase) # BiocManager::install('GSEABase')
  ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
  ## 如果存在之前分析后保存的结果文件，就不需要重复进行GSEA分析。
  f=paste("step3.2_DEG_gsea_results.Rdata")
  if(!file.exists(f)){
    gsea_results <- lapply(gmts, function(gmtfile){
      # gmtfile=gmts[2]
      geneset <- read.gmt(file.path(d,gmtfile)) 
      print(paste0('Now process the ',gmtfile))
      egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
      head(egmt)
      # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
      
      return(egmt)
    })
    # 上面的代码耗时，所以保存结果到本地文件
    save(gsea_results,file = f)
  }
  load(file = f)
  #提取gsea结果，熟悉这个对象
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  write.csv(gsea_results_df,file = paste("step3.2_DEG_gsea_results.csv"))
  
  pdf(paste("step3.2_Top3_pathway_gsea.pdf"),width = 10,useDingbats = F)
  library(RColorBrewer)
  mycol=brewer.pal(4,"Set2")
  gseaplot2(gsea_results[[1]],c("HALLMARK_APOPTOSIS","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),color = mycol,pvalue_table = F) 
  dev.off()
  
  
}
