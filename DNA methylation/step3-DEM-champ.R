rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)

load(file = "step2-champ_myNorm.Rdata")
myLoad    # 存储了甲基化信号矩阵和表型信息。



group_list=myLoad$pd$Sample_Group
table(group_list)
myDMP <- champ.DMP(beta = myNorm2,pheno=group_list,adjPVal = 1,compare.group = c("AMD","Control"))
# champ.DMP函数用于分析差异甲基化探针，champ.DMR函数用于分析差异甲基化区域。
head(myDMP[[1]])
# 测试数据只有两个分组，所以list 中只有一个元素。差异分析的结果是一个data.frame对象，可以分成3个部分。
# 从logFC到B的部分是limma 差异输出结果， C_AVG到deltaBeta是每组表达量的均值，
# deltaBate是两组均值的差，CHR到Probe_SNPs_10是探针的注释信息。
save(myDMP,file = 'step3-output-myDMP.Rdata')




