rm(list = ls())   
options(stringsAsFactors = F)

load(file="step3.1-sub_deg.Rdata")
library(dplyr)
library(plyr)
library(ggplot2)
### 染色体分布图
if(F){
  # BiocManager::install("RIdeogram")
  library(RIdeogram)
  data(human_karyotype, package="RIdeogram")
  data(gene_density, package="RIdeogram")
  data(Random_RNAs_500, package="RIdeogram")
  head(human_karyotype)
  head(gene_density)
  head(Random_RNAs_500)
  
  
  
  
  
  
  
  .library(stringr)
  {
    tmp=str_match(sub_deg$UCSC_CpG_Islands_Name,"chr([0-9]+):([0-9]+)-([0-9]+)")
    gene_density=data.frame(Chr=as.numeric(tmp[,2]),Start=as.numeric(tmp[,3]),End=as.numeric(tmp[,4]),Value=sub_deg$deltaBeta)
    ##有很多没给范围的
    gene_density=data.frame(Chr=as.numeric(sub_deg$CHR),Start=as.numeric(sub_deg$MAPINFO),End=as.numeric(sub_deg$MAPINFO)+1,Value=sub_deg$deltaBeta)
    ideogram(karyotype = human_karyotype, overlaid = gene_density)
    convertSVG("chromosome.svg", device = "png")
    svg2pdf("chromosome.svg")
  }###线太细了 
  
  gene_density=data.frame(Chr=as.numeric(sub_deg$CHR),Start=as.numeric(sub_deg$MAPINFO),End=as.numeric(sub_deg$MAPINFO),Value=sub_deg$betadif)
  
  
  
  {###以下是
    library(chromoMap)
    write.table(human_karyotype[,1:4],file = "chromosome_file.txt",row.names = F,col.names = F,quote = F,sep = "\t")
    annotation_file=data.frame(name=rownames(sub_deg),Chr=as.character(sub_deg$CHR),Start=as.numeric(sub_deg$MAPINFO),End=as.numeric(sub_deg$MAPINFO)+1,Value=sub_deg$deltaBeta)
    
    write.table(annotation_file,file = "annotation_file.txt",row.names = F,col.names = F,quote=F,sep = "\t")

    chromoMap("chromosome_file.txt","annotation_file.txt",
              data_based_color_map = T,
              chr_color = c("black"),
              anno_col = c("red","blue"),
              chr_width = 10,
              chr_length = 4,
              ch_gap = 5,
              data_type = "numeric",
              v_align =T,
              legend = T,
              lg_x = 100,
              lg_y = 100
              # top_margin = 1,
              # left_margin = 10
    )
    
  }
}



