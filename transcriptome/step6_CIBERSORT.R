rm(list = ls())
options(stringsAsFactors = F)
load("D:/bioinformation/R/甲基化/AMD/转录组/step1_count_anno.Rdata")


library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)
gene_map<-mapIds(org.Hs.eg.db, as.character(rownames(counts)), 'SYMBOL', keytype ='ENSEMBL',multiVals = "first")
gene_map=data.frame(ENSEMBL=names(gene_map),SYMBOL=as.character(gene_map))



# 因为CIBERSORT推荐是用TPM做输入文件，所以需要将count 转换为TPM
{
  # 载入GTF文件
  library(GenomicFeatures)
  txdb <- makeTxDbFromGFF("D:/bioinformation/GTF注释文件/Homo_sapiens.GRCh38.94.gtf",format="gtf")
  # 通过exonsBy获取每个gene上的所有外显子的起始位点和终止位点，然后用reduce去除掉重叠冗余的部分，最后计算长度
  
  exons_gene <- exonsBy(txdb, by = "gene")
  exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
  exons_gene_lens=t(as.data.frame(exons_gene_lens))
  counts_length=merge(counts,exons_gene_lens,by="row.names")
  rownames(counts_length)=counts_length[,1]
  counts_length=counts_length[,-1]
  # TPM计算
  
  kb <- counts_length$V1 / 1000
  kb
  countdata <- counts_length[,-ncol(counts_length)]
  rpk <- countdata / kb
  rpk
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  tpm=as.data.frame(tpm)
  head(tpm)
  save(tpm,file = "step1_tpm.data")
  # write.table(tpm,file="step1_tpm.xls",sep="\t",quote=F)
  # FPKM计算
  fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
  head(fpkm)
  write.table(fpkm,file="fpkm.xls",sep="\t",quote=F)
  
  # FPKM转化为TPM
  fpkm_to_tpm = t(t(fpkm)/colSums(fpkm))*10^6
  head(fpkm_to_tpm)
  
  
 
}

dat =tpm
dat$median=apply(dat,1,median)
dat$ENSEMBL = rownames(dat)
dat = merge(dat,gene_map,by="ENSEMBL")
dat = na.omit(dat)

dat=dat[order(dat$SYMBOL,dat$median,decreasing = T),]#把dat$symbol按照dat$median排序
dat=dat[!duplicated(dat$SYMBOL),]#取出不重复的dat$symbol
rownames(dat)=dat$SYMBOL
dat=dat[,-c(1,ncol(dat)-1,ncol(dat))]
dim(dat) 
write.table(dat,file = "step1_TPM_SYMBOL.txt",quote = F,sep="\t")
# 第一行第一列交界处需要用Tab键隔开一格，不然会导致报错欧。

## CIBERSORT
{
  source("D:/bioinformation/R/各种功能/CIBERSORT/CIBERSORT.R")
  sig_matrix <- "D:/bioinformation/R/各种功能/CIBERSORT/LM22.txt"
  ## 指定表达矩阵
  mixture_file = 'step1_TPM_SYMBOL.txt'
  # 运行
  res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=10, QN=F)
  # 其中nperm给的是置换的次数，QN如果是芯片设置为T，如果是测序就设置为F，测序数据最好是TPM
  # 文件夹中会出现CIBERSORT-Results.txt文件，即为结果
  
  #可视化
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  cibersort_raw <- cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
    rename("Mixture" = "Sample") %>%
    dplyr::select(-c("P.value","Correlation","RMSE"))
  # 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Sample”。
  #并赋值给cibersort_raw。
  
  cibersort_tidy <- cibersort_raw %>%
    remove_rownames() %>%
    column_to_rownames("Sample")
  # 将cibersort_raw第一列变为列名后赋值给cibersort_tidy。
  
  flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) < 
                  dim(cibersort_tidy)[1]/2)
  # 筛选出0值太多的一些细胞
  
  cibersort_tidy <- cibersort_raw[,which(flag)] %>%
    as.matrix() %>%
    t()
  # 留下在大部分样本中有所表达的细胞。
  
  bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
  # breaks用来定义数值和颜色的对应关系。
  
  # 将CIBERSORT_Result进行可视化
  #柱状图可视化细胞占比预测
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  cibersort_barplot <- cibersort_raw %>%
    gather(key = Cell_type,value = Proportion,2:23)
  #使用RColorBrewer包配置需要的色彩方案，使用gather函数中的key-value对应关系重建细胞名称和比例的对应关系并赋值给cibersort_barplot
  
  ggplot(cibersort_barplot,aes(Sample,Proportion,fill = Cell_type)) + 
    geom_bar(position = "stack",stat = "identity") +
    labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + theme_bw() +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_y_continuous(expand = c(0.01,0)) +
    scale_fill_manual(values = mypalette(23))+ggsave("step6_堆积图.pdf",width=10,useDingbats = F)
  # 百分比堆积图
  
  
  # 箱线图
  ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
    geom_boxplot(outlier.shape = 21,color = "black") + theme_bw() + 
    labs(x = "Cell_Type", y = "Estimated Proportion") +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = mypalette(23))+ggsave("step6_all_boxplot.pdf",width=10,useDingbats = F)
  #调整参数让柱状图更加美观。
  
  # 按分组箱线图
  require(ggplot2)
  require(ggsci)
  require(ggpubr)
  require(ggthemes)
  plot.info=merge(cibersort_barplot,pdata,by.x = "Sample", by.y = "geo_accession")
  ggboxplot(
      plot.info,
      x = "Cell_type",
      y = "Proportion",
      color = "black",
      fill = "group",
      xlab = "",
      ylab = "Cell composition",
      main = "TME Cell composition group by Group"
    ) +
      stat_compare_means(
          label = "p.signif",
          method = "wilcox.test",
          ref.group = ".all.",
          hide.ns = F
        ) +
      theme_base() +
      theme(axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 1
        ))+ggsave('boxplot.immunecell.pdf',width=16,useDingbats = F)
  ##热图
  {
    cibersort_raw <- read.table("CIBERSORT-Results.txt",header = T,sep = '\t') %>%
      rename("Patients" = "Mixture") %>%
      select(-c("P.value","Correlation","RMSE"))
    cibersort_tidy <- cibersort_raw %>%
      remove_rownames() %>%
      column_to_rownames("Patients")
    
    
    flag <- apply(cibersort_tidy,2,function(x) sum(x == 0) < 
                    dim(cibersort_tidy)[1]/2)
    
    
    cibersort_tidy <- cibersort_tidy[,which(flag)] %>%
      as.matrix() %>%
      t()
    
    
    bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
    
    
    
    
    library(pheatmap)
    library(RColorBrewer)
    pheatmap(
      cibersort_tidy,
      breaks = bk,
      cluster_cols = T,
      scale = "row",
      cluster_row = T,
      border_color = NA,
      show_colnames = F,
      show_rownames = T,
      color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                colorRampPalette(colors = c("white","red"))(length(bk)/2)
      ))
    
    }
}

