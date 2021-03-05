rm(list = ls())
options(stringsAsFactors = F)

library(ArrayExpress)
AE = getAE("E-MTAB-7183")
save(AE,file = "AE.Rdata")
pdata=read.table("E-MTAB-7183.sdrf.txt",sep = "\t",header = T)##临床信息

library("ChAMP")#450测试数据testDir=system.file("extdata",package="ChAMPdata")
dir= "D:/bioinformation/R/甲基化/AMD/甲基化/Idat"
myimport <- champ.import(dir,arraytype="450K")#
myimport
dim(myimport$beta)
myLoad=champ.filter(beta=myimport$beta,pd=myimport$pd)
dim(myLoad$beta)
Sample_Group=myLoad$pd$Sample_Group
Sample_Group=ifelse(Sample_Group=="normal","Control","AMD")
myLoad$pd$Sample_Group=Sample_Group
champ.QC(beta=myLoad$beta,pheno=myLoad$pd$Sample_Group)
QC.GUI()
save(myimport,myLoad,file = 'step1-output.Rdata')


