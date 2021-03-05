rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')

myLoad  
dim(myLoad$beta)
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=3)
dim(myNorm) 
champ.SVD(beta = myNorm,pd=myLoad$pd)
champ.QC(beta=myNorm,pheno=myLoad$pd$Sample_Group,resultsDir="./CHAMP_NormQCimages/")
pdata=myLoad$pd
# MDSplot
{
  library(limma)
  # 筛选top1000 var
  rv=apply( myNorm, 1, var)
  expr=head(myNorm[order(rv,decreasing = T),],1000)
  
    # Classical MDS
  # N rows (objects) x p columns (variables)
  # each row identified by a unique row name
  
  d <- dist(t(expr)) # euclidean distances between the rows
  fit <- cmdscale(d, eig=T,k=2) # k is the number of dim
  fit # view results
  # plot solution 
  data=cbind(x=fit$points[,1],y=fit$points[,2],pdata)
  ggplot(data,aes(x=x,y=y,colour=Sample_Group))+
    geom_point(alpha=0.5)+
    labs(x = "Coordinate 1",y="Coordinate 2")+
    ggtitle("Top 1000 most variable Cpgs") +
    scale_colour_brewer(palette = "Set1")+
    theme(plot.title = element_text(hjust = 0.5))+##标题居中
  ggsave("Top 1000 most variable Cpgs-Group.pdf",width =4,height = 4, useDingbats=FALSE) 
  
  ggplot(data,aes(x=x,y=y,colour=Slide))+
    geom_point(alpha=0.5)+
    labs(x = "Coordinate 1",y="Coordinate 2")+
    ggtitle("Top 1000 most variable Cpgs") +
    scale_colour_brewer(palette = "Set1")+
    theme(plot.title = element_text(hjust = 0.5))+##标题居中
  ggsave("Top 1000 most variable Cpgs-Slide.pdf",width =4,height = 4, useDingbats=FALSE) 
  
  ggplot(data,aes(x=x,y=y,colour=sex))+
    geom_point(alpha=0.5)+
    labs(x = "Coordinate 1",y="Coordinate 2")+
    ggtitle("Top 1000 most variable Cpgs") +
    scale_colour_brewer(palette = "Set1")+
    theme(plot.title = element_text(hjust = 0.5))+##标题居中
    ggsave("Top 1000 most variable Cpgs-Sex.pdf",width =4,height = 4, useDingbats=FALSE) 
  
  plot(x, y, xlab="", ylab="Coordinate 2", 
       main="Metric MDS", type="n")
  text(x, y, labels = row.names(mydata), cex=.7)  
}
# 变异系数解释度
if(F){
  
  library('variancePartition')
  info=myLoad$pd
  geneExpr=myLoad$beta
  rv=apply( geneExpr, 1, var)
  expr=head(myNorm[order(rv,decreasing = T),],1000)
  geneExpr=expr
  str(info)
  form <- ~ age + (1|sex) + (1|Slide) + (1|Array) + (1|Sample_Group)
  varPart <- fitExtractVarPartModel( geneExpr, form, info )
  # sort variables (i.e. columns) by median fraction
  # of variance explained
  vp <- sortCols( varPart )
  
  # Figure 1a
  # Bar plot of variance fractions for the first 10 genes
  plotPercentBars( vp[1:10,] )
  ggsave("Top 10 variance explained.pdf",width =4,height = 4, useDingbats=FALSE) 
  #
  # Figure 1b
  # violin plot of contribution of each variable to total variance
  plotVarPart( vp )
  ggsave("variance explained-TOP1000.pdf",width =4,height = 4, useDingbats=FALSE)
  # Access first entries
  head(varPart)
  # Access first entries for Individual
  head(varPart$Slide)
  # sort genes based on variance explained by Individual
  head(varPart[order(varPart$Individual, decreasing=TRUE),])
}


str(myLoad$pd)
myLoad$pd$Slide=as.factor(myLoad$pd$Slide)


myNorm2=champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname = c("Slide"))
champ.SVD(beta = myNorm2,pd=myLoad$pd,resultsDir="./CHAMP_SVDimages_combat/")
champ.QC(beta=myNorm2,pheno=myLoad$pd$Sample_Group,resultsDir="./CHAMP_QCimages__combat/")
save(myNorm,myLoad,myNorm2,file = 'step2-champ_myNorm.Rdata')


