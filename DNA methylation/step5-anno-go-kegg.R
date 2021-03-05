
rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library(stringr)
library("minfi")
require(GEOquery)
require(Biobase)

load(file = 'step3-output-myDMP.Rdata')
deg=myDMP[[1]]
head(deg)
length(unique(deg$gene))

deltabeta_cutoff=0
deg$change = as.factor(ifelse(deg$P.Value < 0.01 & abs(deg$deltaBeta) > deltabeta_cutoff,
                                   ifelse(deg$deltaBeta > deltabeta_cutoff ,'UP','DOWN'),'NOT'))


table(deg$change)
head(deg)
sub_deg = deg[deg$change %in% c("UP","DOWN"),]
write.csv(sub_deg,file = "sub_deg.csv")
sub_deg$symbol=sub_deg$gene
sub_deg=sub_deg[!sub_deg$feature=="IGR",]
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(sub_deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)



head(sub_deg)

sub_deg=merge(sub_deg,df,by.y='SYMBOL',by.x='symbol')
head(sub_deg)
save(sub_deg,file = 'anno_DEG.Rdata')


gene_up= sub_deg[sub_deg$change == 'UP','ENTREZID'] 
gene_down=sub_deg[sub_deg$change == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)




gene_up=unique(gene_up)
gene_down=unique(gene_down)
gene_diff=unique(c(gene_up,gene_down))
pro='test_methy'
###   over-representation test
# 下面把3个基因集分开做超几何分布检验
# 首先是上调基因集。
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
kk=kk.up
dotplot(kk)
tmp=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv(tmp@result,'GO_diff_all.csv')
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv(kk@result,paste0(pro,'_kk.up.csv'))

# 首先是下调基因集。
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
head(kk.down)[,1:6]
kk=kk.down
dotplot(kk)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv(kk@result,paste0(pro,'_kk.down.csv'))

# 最后是上下调合并后的基因集。
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'hsa',
                      pvalueCutoff = 0.9)
head(kk.diff)[,1:6]
kk=kk.diff
dotplot(kk, font.size =10)+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave(filename = paste0(pro,'_kegg_up_down.pdf'),height = 6,width=8,useDingbats = F)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv(kk@result,paste0(pro,'_kk.diff.csv'))


kegg_diff_dt <- as.data.frame(kk.diff)
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)

ggsave(g_kegg,filename = paste0(pro,'_kegg_up_down.png') )



go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_up_GO_all_barplot.pdf', useDingbats=FALSE) 
  tmp=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(tmp@result,'GO_up_all.csv')
  

go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_down_GO_all_barplot.pdf', useDingbats=FALSE)
  tmp=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(tmp@result,'GO_down_all.csv')
  
  
go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
  ggsave('gene_diff_GO_all_barplot.pdf', useDingbats=FALSE)
  tmp=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(tmp@result,'GO_diff_all.csv')
  