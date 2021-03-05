rm(list = ls()) 
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

### 对 MigDB中的全部基因集 做GSVA分析。
## 还有ssGSEA, PGSEA
{
  gsea_results_df= read.csv("step3.2_DEG_gsea_results.csv",row.names = 1)
  setnames=rownames(gsea_results_df)
  load(file = "step2-DEG-macular.Rdata")
  # 每次都要检测数据
  group_list=pd$group
  exprSet[1:4,1:4]  
  #ENSEMBL-symbol
  {
    library(org.Hs.eg.db)
    gene_map<-mapIds(org.Hs.eg.db, as.character(rownames(exprSet)), 'SYMBOL', keytype ='ENSEMBL',multiVals = "first")
    exprSet$SYMBOL=gene_map
    exprSet=na.omit(exprSet)
    exprSet=exprSet[!duplicated(exprSet$SYMBOL),]
    exprSet <- exprSet %>% remove_rownames() %>% column_to_rownames("SYMBOL")
    exprSet[1:4,1:4]
  }
  
  X=as.matrix(exprSet)
  
  table(group_list)
  ## Molecular Signatures Database (MSigDb) 
  d='D:/bioinformation/MsigDB/symbols/'
  gmts=list.files(d,pattern = 'all')
  gmts
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSEABase)
  library(GSVA) # BiocManager::install('GSVA')

  
   gmts=c("c2.cp.kegg.v7.1.symbols.gmt" ,"h.all.v7.1.symbols.gmt")
  
  if(T){
    # GSVA将表达矩阵转换成通路富集分数(ES)矩阵
    es_max_all <- lapply(gmts, function(gmtfile){ 
      #gmtfile=gmts[8];gmtfile
      geneset <- getGmt(file.path(d,gmtfile))  
      es.max <- GSVA::gsva(X, geneset, 
                           mx.diff=FALSE, verbose=FALSE, 
                           parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max_all, function(es.max){
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                                     levels = design)
      contrast.matrix<-makeContrasts("AMD-Control",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
  } 
  
  gmts
  
  save(es_max_all,es_deg,file='step3.3_gsva_msigdb.Rdata')
  
  
  load(file='step3.3_gsva_msigdb.Rdata')
  
  
  dat=es_max_all[[2]]
  gsea_gsva_result=dat[setnames,]
  write.csv(gsea_gsva_result,file = 'step3.3-gsea_gsva_result.csv')
  }
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max_all)
  df=df[df$P.Value<0.05 & abs(df$logFC) > 0.1,]
  write.csv(df,file = 'GSVA_DEG.csv')
  
  
  # GSVA&CIBERSORT
  
  {
    library(dplyr)
    library(tidyr)
    library(tidyverse)
    library(ggstatsplot)
    cibersort_raw <- cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
      rename("Mixture" = "Sample") %>%
      dplyr::select(-c("P.value","Correlation","RMSE"))
    # 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Sample”。
    #并赋值给cibersort_raw。
    
    cibersort_tidy <- cibersort_raw %>%
      remove_rownames() %>%
      column_to_rownames("Sample")
    
    
      cibersort_tidy <- cibersort_tidy[colnames(gsea_gsva_result),]
      
      
    
    ##TNFA_SIGNALING_VIA_NFKB
    {
      index <-"HALLMARK_TNFA_SIGNALING_VIA_NFKB" #基因名
      index_expr=gsea_gsva_result[index,]
      y <- as.numeric(index_expr)
      head(y)
      
      colnames <- colnames(cibersort_tidy)
      data <- data.frame(colnames)
      for (i in 1:length(colnames)){
        test <- cor.test(as.numeric(cibersort_tidy[,i]),y, method="spearman")
        data[i,2] <- test$estimate                                            
        data[i,3] <- test$p.value
      }
      names(data) <- c("symbol","correlation","pvalue")
      head(data)
      
      #输出到文件
      write.table(data, "step3.3_TNFA_SIGNALING_VIA_NFKB_cor.txt", sep = "\t", quote = F, row.names = F)
      
      
      head(data)
      data=na.omit(data)
      
      data %>% 
        #filter(pvalue <0.05) %>% # 如果不想把p值大于0.05的放在图上，去掉最前面的#号
        ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
        geom_segment(aes(xend=0,yend=symbol)) +
        geom_point(aes(col=pvalue,size=abs(correlation))) +
        scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
        # scale_color_viridis_c(begin = 0.5, end = 1) +
        scale_size_continuous(range =c(2,8))  +
        theme_minimal() +
        ylab(NULL)+
      
      ggsave("step3.3_TNFA_SIGNALING_VIA_NFKB_Xcell.pdf")
      
      
      
      # 挑选细胞跟表达量散点图
      imucell <- "Macrophages.M1"
      # 合并免疫数据和表达量数据
      plot_df <- data.frame(gene=y,imucell=cibersort_tidy[,imucell])
      head(plot_df)
      
      pdf("step3.3_TNFA_SIGNALING_VIA_NFKB_X_M1.pdf",useDingbats = F)
      ggscatterstats(data = plot_df,
                     x = gene,
                     y = imucell,
                     centrality.para = NULL,
                     margins = "both",
                     
                     xfill = "#CC79A7",
                     yfill = "#009E73"
                     ,marginal.type = NULL
                     )
      dev.off()
      
    }
      
      ##HALLMARK_INFLAMMATORY_RESPONSE"
      {
        index <-"HALLMARK_INFLAMMATORY_RESPONSE" #基因名
        index_expr=gsea_gsva_result[index,]
        y <- as.numeric(index_expr)
        head(y)
        
        colnames <- colnames(cibersort_tidy)
        data <- data.frame(colnames)
        for (i in 1:length(colnames)){
          test <- cor.test(as.numeric(cibersort_tidy[,i]),y, method="spearman")
          data[i,2] <- test$estimate                                            
          data[i,3] <- test$p.value
        }
        names(data) <- c("symbol","correlation","pvalue")
        head(data)
        
        #输出到文件
        write.table(data, "step7_INFLAMMATORY_RESPONSE_cor.txt", sep = "\t", quote = F, row.names = F)
        
        
        head(data)
        data=na.omit(data)
        
        data %>% 
          #filter(pvalue <0.05) %>% # 如果不想把p值大于0.05的放在图上，去掉最前面的#号
          ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
          geom_segment(aes(xend=0,yend=symbol)) +
          geom_point(aes(col=pvalue,size=abs(correlation))) +
          scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
          # scale_color_viridis_c(begin = 0.5, end = 1) +
          scale_size_continuous(range =c(2,8))  +
          theme_minimal() +
          ylab(NULL)+
          
          ggsave("step7_INFLAMMATORY_RESPONSE_Xcell.pdf")
      }
  }
  
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=1
    print(i)
    dat=es_max_all[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.05,]
    print(dim(df))
    if(nrow(df)>5){
      n=rownames(df)
      dat=dat[match(n,rownames(dat)),]
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(dat)
      rownames(dat)=substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
      
    }
  })
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max_all)
  df=df[df$P.Value<0.05 & abs(df$logFC) > 0.1,]
  write.csv(df,file = 'GSVA_DEG.csv')
}





