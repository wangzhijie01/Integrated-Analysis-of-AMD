rm(list = ls())
options(stringsAsFactors = F)

load(file="step3-output-myDMP.Rdata")
library(dplyr)
library(plyr)
library(ggplot2)

{
  deg=myDMP[[1]]
  head(deg)
  length(unique(deg$gene))
  
  deltabeta_cutoff=0
  deg$change = as.factor(ifelse(deg$P.Value < 0.01 & abs(deg$deltaBeta) > deltabeta_cutoff,
                                ifelse(deg$deltaBeta > deltabeta_cutoff ,'UP','DOWN'),'NOT'))
  
  
  table(deg$change)
  head(deg)
  sub_deg = deg[deg$change %in% c("UP","DOWN"),]
  sub_deg$cpg=rownames(sub_deg)
    dim(sub_deg)
  library(ggplot2)
  save(sub_deg,file="step3.1-sub_deg.Rdata")
  
}





library(dplyr)
data_group<- group_by(sub_deg, change,feature)
data_GroupByID<- dplyr::summarise(data_group,count = n())
data_GroupByID = plyr::ddply(data_GroupByID,"feature",transform,label_y=cumsum(count)-0.5*count)
# data_GroupByID<- data_GroupByID[order(data_GroupByID$count,decreasing=T),]#降序排序


{
  
  ggplot(data_GroupByID, aes(x = feature, y = count, fill = change, label = count)) +
    geom_bar(stat = "identity",colour="black") +
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_manual(values=c( "#b3cde3","#fbb4ae"), 
                      breaks=c("DOWN", "UP"),
                      labels=c( "hypomethylated","hypermethylated"))
  ggsave("DMP-feature-柱形图.pdf",width = 8,height = 4, useDingbats=FALSE)
  
}


{
  data_Group<- dplyr::group_by(sub_deg, change,cgi)
  data_GroupByID<- dplyr::summarise(data_Group,count = n())
  data_GroupByID= plyr::ddply(data_GroupByID,"cgi",transform,label_y=cumsum(count)-0.5*count)
  # data_GroupByID<- data_GroupByID[order(data_GroupByID$count,decreasing=T),]#降序排序
  ggplot(data_GroupByID,aes(x=cgi,y=count,fill=change, label = count))+
    geom_bar(stat="identity",colour="black")+
    geom_text(size = 3, position = position_stack(vjust = 0.5))+
    scale_fill_manual(values=c( "#b3cde3","#fbb4ae"), 
                      breaks=c("DOWN", "UP"),
                     labels=c( "hypomethylated","hypermethylated"))
    # scale_fill_brewer(palette = "Pastel1")
  ggsave("DMP-cgi-柱形图.pdf",width = 6,height = 4, useDingbats=FALSE)
}

##注释
{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df <- bitr(unique(sub_deg$gene), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  sub_deg=merge(sub_deg,df,by.x="gene",by.y="SYMBOL")
  dim(sub_deg)
  table(sub_deg$change)
 
}
