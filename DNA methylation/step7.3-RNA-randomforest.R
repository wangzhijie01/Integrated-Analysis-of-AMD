
rm(list = ls())   
options(stringsAsFactors = F)

load("step7-DEG.DMG.Rdata")##merge
load("D:/bioinformation/R/甲基化/AMD/转录组/step1_count_anno.Rdata")##转录组表达量
RNA.pdata=pdata

###数据预处理###

if (F) {

  gene=unique(c(as.character(hypo_up$gene),as.character(hyper_down$gene)))
  ENSEMBL=unique(c(as.character(hypo_up$ENSEMBL),as.character(hyper_down$ENSEMBL)))
  ids=data.frame(gene=gene,ENSEMBL=ENSEMBL)
  # 提取macular的数据
  {
    pd=RNA.pdata[grep("RPE, Macula",RNA.pdata$tissue),]
    dim(pd)
    table(pd$tissue,pd$group)
    exprSet=counts[,pd$geo_accession]
    dim(exprSet)
  }
  
  RNAseq=as.data.frame(t(exprSet[ENSEMBL,]))
  colnames(RNAseq)=gene
  RNAseq=log2(RNAseq+1)
  dim(RNAseq)
  RNAseq$disease=as.factor(pd$group)
  dim(RNAseq)
  str(RNAseq)
}


##randome-forest 建模##

library(randomForest)


##参数筛选##
if(T){
  # 筛选ntree
  set.seed(100)
  rf_ntree <- randomForest(disease~.,data=RNAseq,ntree=1000)
  pdf("RNA_ntree_chosen.pdf",useDingbats = F)
  plot(rf_ntree)
  # 200
  dev.off()
  # 筛选mtry
  n <- ncol(RNAseq)-1
  errRate <- c()
  
  # 从图中可以看到，当ntree=400时，模型内的误差就基本稳定了，出于更保险的考虑，我们确定ntree值为100。
  set.seed(100)
  for (i in 1:n){ 
    
    m <- randomForest(disease~.,data=RNAseq,mtry=i,proximity=TRUE) 
    
    err<-mean(m$err.rate[,1])
    
    errRate[i] <- err
  }  
  print(errRate)
  m= which.min(errRate)  
  print(m)
  # [1] 10
  pdf("RNA_ntry_chosen.pdf",useDingbats = F)
  matplot(1:n , errRate, pch=19 , col=c("red"),type="p",ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
  legend("topright",legend=c("Out of Bag Error"),pch=19, col=c("red"))
  dev.off()
  #选择最优mtry参数值
  
}

##importance 筛选基因排序##
if (F) {
  set.seed(100)
  rf_output=randomForest(disease~., data =RNAseq ,importance = TRUE,mtry=10,ntree=200,
                         proximity=TRUE)
  rf_importances=importance(rf_output, scale=FALSE)
  choose_gene=rownames(tail(rf_importances[,],10))
  head(rf_importances)
  gene.order=rownames(rf_importances)[order(rf_importances[,4],decreasing = T)]
  write.csv(rf_importances,file = "RNA_importance.csv")
  
  pdf("RNAseq-rf_importance.pdf",useDingbats = F)
  varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, 
             main="Variable Importance (Gini) for top 30 predictors",cex =.7)
  dev.off()
}


##caret包train函数测试##
if (F) {
  # Leave One Out Cross Validation
  # load the library
  library(caret)
  # define training control
  train_control <- trainControl(method="LOOCV",
                                summaryFunction=twoClassSummary, classProbs=TRUE,##把评价改成ROC，不然是accuracy
                                savePredictions = T)
  grid=expand.grid(mtry=10)
  # train the model
  model_RNA <- train(disease~., data =RNAseq,trControl=train_control,method="rf",tuneGrid =grid,ntree=200)
  # summarize results
  print(model_RNA)
  model_RNA$results
}


##按importance顺序逐个纳入基因，寻找最佳纳入基因数量##
# set.seed(12345679)
# sam<- createDataPartition(methylation$disease, p = .5,list = FALSE)
# train <- methylation[sam,]
# test <- methylation[-sam,]
if (T) {
  ROC=c()
  RNA_model=list()
  for (i in 1:length(gene.order)) {
    set.seed(12345679)
    gen.enroll = gene.order[1:i]
    tmp.data = RNAseq[,c(gen.enroll,"disease")]
    #Leave One Out Cross Validation
    train_control <- trainControl(method="LOOCV",
                                  summaryFunction=twoClassSummary, classProbs=TRUE,##把评价改成ROC，不然是accuracy
                                  savePredictions = T)
    
    #Repeated k-fold Cross Validation
    # train_control <- trainControl(method="repeatedcv", number=10, repeats=3,
    #                               summaryFunction=twoClassSummary, classProbs=TRUE,##把评价改成ROC，不然是accuracy
    #                               savePredictions = T)##把predict值存储下来画ROC用
    grid = expand.grid(mtry=10)
    # train the model
    model_RNA <- train(disease~., data =tmp.data,trControl=train_control,method="rf",tuneGrid =grid,ntree=200)
    # summarize results
    print(model_RNA)
    ROC[i] = model_RNA$results$ROC[1]
    RNA_model[[i]] = model_RNA
  }
  
  
  ROC
  
  pdf("Gene expresion classifier_chosen_gene.pdf",useDingbats = F)
  plot(1:n,ROC, pch=19 , col=c("black"),type="p",cex=1,
       ylab="AUC of ROC curve",xlab="Number of Predictors in classifier",
       main="Gene expresion classifier")
  dev.off()

  ntmp=which.max(ROC)
  ntmp
  final_model = RNA_model[[ntmp]]
  ROC[ntmp]
  save(ROC,RNA_model,final_model,ids,file = "step7.3-RNAseq-final_rf_model.Rdata")
  
  library(pROC)
  pred.tmp=model_RNA$pred$AMD
  y.tmp=model_RNA$pred$obs
  plot.roc(y.tmp,
           pred.tmp)
  #P值和AUC
  {
    rf.ROC = roc(response = y.tmp,
                 predictor = pred.tmp)
    auc=auc(rf.ROC)
    library(verification)
    y1=ifelse(y.tmp=="AMD",1,0)
    ROC.test=roc.area(y1 ,pred.tmp)
    P.value=ROC.test$p.value
  }
  
  
  ##画ROC曲线
  pdf("RNA-ROC_train.pdf",useDingbats = F)
  title=paste("Top",ntmp,"gene expression classifier")
  plot(rf.ROC,
       print.auc=TRUE, print.auc.x=0.4, print.auc.y=0.5,
       # 图像上输出AUC值,坐标为（x，y）
       auc.polygon=TRUE, auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色
       max.auc.polygon=FALSE,  # 填充整个图像
       grid=c(0.5, 0.2), grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色
       print.thres=TRUE, print.thres.cex=0.9, # 图像上输出最佳截断值，字体缩放倍数
       smooth=F, # 绘制不平滑曲线
       main=title, # 添加标题
       col="#FF2E63",  # 曲线颜色
       legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度
  text(0.4, 0.4, labels=paste("P value =", format.pval(P.value)), col="#FF2E63",adj=c(0, .5)) # 在图上添加P值
  dev.off()
}


# 5.模型预测和评估


load("step7.3-RNAseq-final_rf_model.Rdata")
# 1.首先用non macular 的数据验证
# 提取macular的数据
{
  load("D:/bioinformation/R/甲基化/AMD/转录组/step1_count_anno.Rdata")##转录组表达量
  RNA.pdata=pdata
  {
    pd=RNA.pdata[grep("RPE, non-Macula",RNA.pdata$tissue),]
    dim(pd)
    table(pd$tissue,pd$group)
    exprSet.non=counts[,pd$geo_accession]
    dim(exprSet.non)
  }
  
  RNAseq.non=as.data.frame(t(exprSet.non[ENSEMBL,]))
  colnames(RNAseq.non)=gene
  RNAseq.non=log2(RNAseq.non+1)
  dim(RNAseq.non)
  RNAseq.non$disease=as.factor(pd$group)
  dim(RNAseq.non)
  str(RNAseq.non)
  
  test=RNAseq.non
  
  y=RNAseq.non$disease
  rf.prob <- predict(final_model, test,type = "prob")
  
  
  rf.ROC = roc(response = y,
                predictor = rf.prob$AMD)
  auc=auc(rf.ROC)
  ##为计算P值
  {
    library(verification)
    y1=ifelse(y=="AMD",1,0)
    ROC.test=roc.area(y1 ,rf.prob$AMD)
    P.value=ROC.test$p.value
  }
  pdf("RNA-ROC_non_macular.pdf",useDingbats = F)
  title=paste("Top",ntmp,"gene expression classifier")
  plot(rf.ROC,
       print.auc=TRUE, print.auc.x=0.4, print.auc.y=0.5,
       # 图像上输出AUC值,坐标为（x，y）
       auc.polygon=TRUE, auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色
       max.auc.polygon=FALSE,  # 填充整个图像
       grid=c(0.5, 0.2), grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色
       print.thres=TRUE, print.thres.cex=0.9, # 图像上输出最佳截断值，字体缩放倍数
       smooth=F, # 绘制不平滑曲线
       main=title, # 添加标题
       col="#FF2E63",  # 曲线颜色
       legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度
  text(0.4, 0.4, labels=paste("P value =", format.pval(P.value)),col="#FF2E63", adj=c(0, .5)) # 在图上添加P值
  dev.off()
  result.coords <- coords(rf.ROC, "best", best.method="closest.topleft", ret=c("all"))
  result.coords
  # plot(rf.ROC,type = "S",col = "red")
  # text(0.8,0.2, labels = paste0("AUC = ",round(auc,3)))#> [1] 1
}

rm(list = ls())
# 2.用macular RPE转录组芯片进行验证
load("step7.3-RNAseq-final_rf_model.Rdata")
load('D:/bioinformation/R/甲基化/AMD/验证数据集/GSE29801/step1_count_anno.Rdata')

{
  
  counts[1:4,1:4]
  RNAseq.non=as.data.frame(t(counts[gene,]))
  RNAseq.non=log2(RNAseq.non+1)
  dim(RNAseq.non)
  RNAseq.non$disease=as.factor(pdata$group)
  dim(RNAseq.non)
  str(RNAseq.non)
  
  test=RNAseq.non
  
  y=RNAseq.non$disease
  rf.prob <- predict(final_model, test,type = "prob")
  
  
  rf.ROC = roc(response = y,
               predictor = rf.prob$AMD)
  auc=auc(rf.ROC)
  ##为计算P值
  {
    library(verification)
    y1=ifelse(y=="AMD",1,0)
    ROC.test=roc.area(y1 ,rf.prob$AMD)
    P.value=ROC.test$p.value
  }
  pdf("RNA-ROC_GSE29801.pdf",useDingbats = F)
  title=paste("gene expression classifier")
  plot(rf.ROC,
       print.auc=TRUE, print.auc.x=0.4, print.auc.y=0.5,
       # 图像上输出AUC值,坐标为（x，y）
       auc.polygon=TRUE, auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色
       max.auc.polygon=FALSE,  # 填充整个图像
       grid=c(0.5, 0.2), grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色
       print.thres=TRUE, print.thres.cex=0.9, # 图像上输出最佳截断值，字体缩放倍数
       smooth=F, # 绘制不平滑曲线
       main=title, # 添加标题
       col="#FF2E63",  # 曲线颜色
       legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度
  text(0.4, 0.4, labels=paste("P value =", format.pval(P.value)), adj=c(0, .5)) # 在图上添加P值
  dev.off()
  result.coords <- coords(rf.ROC, "best", best.method="closest.topleft", ret=c("all"))
  result.coords
  # plot(rf.ROC,type = "S",col = "red")
  # text(0.8,0.2, labels = paste0("AUC = ",round(auc,3)))#> [1] 1
}





# 3.用500+retina的测序转录组数据
rm(list = ls())
load("step7.3-RNAseq-final_rf_model.Rdata")

if (F) {
  load("D:/bioinformation/R/甲基化/AMD/验证数据集/GSE115828/step1_ensembl.Rdata")###AMD转录组验证集,已取log2
  pdata = pdata[!pdata$mgs_level==2,]
  counts_expr=counts_expr[,rownames(pdata)]
  counts_expr[1:4,1:4]
  ENSEMBL=ids$ENSEMBL
  RNAseq.non=as.data.frame(t(counts_expr[ENSEMBL,]))
  colnames(RNAseq.non)=ids$gene
  RNAseq.non=log2(RNAseq.non+1)
  dim(RNAseq.non)
  RNAseq.non$disease=as.factor(pdata$disease)
  dim(RNAseq.non)
  str(RNAseq.non)
  
  test=RNAseq.non
  
  y=pdata$disease
  rf.prob <- predict(final_model, test,type = "prob")
  
  
  rf.ROC = roc(response = y,
               predictor = rf.prob$AMD)
  auc=auc(rf.ROC)
  ##为计算P值
  {
    library(verification)
    y1=ifelse(y=="AMD",1,0)
    ROC.test=roc.area(y1 ,rf.prob$AMD)
    P.value=ROC.test$p.value
  }
  pdf("RNA-ROC_GSE115828.pdf",useDingbats = F)
  title=paste("GSE115828 retina samples ROC")
  plot(rf.ROC,
       print.auc=TRUE, print.auc.x=0.4, print.auc.y=0.5,
       # 图像上输出AUC值,坐标为（x，y）
       auc.polygon=TRUE, auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色
       max.auc.polygon=FALSE,  # 填充整个图像
       grid=c(0.5, 0.2), grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色
       print.thres=TRUE, print.thres.cex=0.9, # 图像上输出最佳截断值，字体缩放倍数
       smooth=F, # 绘制不平滑曲线
       main=title, # 添加标题
       col="#FF2E63",  # 曲线颜色
       legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度
  text(0.4, 0.4, labels=paste("P value =", format.pval(P.value)), col="#FF2E63",adj=c(0, .5)) # 在图上添加P值
  dev.off()
  result.coords <- coords(rf.ROC, "best", best.method="closest.topleft", ret=c("all"))
  result.coords
  # plot(rf.ROC,type = "S",col = "red")
  # text(0.8,0.2, labels = paste0("AUC = ",round(auc,3)))#> [1] 1
  
}

# 3.用剩余的retina的测序转录组数据
rm(list = ls())
load("step7.3-RNAseq-final_rf_model.Rdata")
if (F) {
  load("D:/bioinformation/R/甲基化/AMD/转录组/step1-retina.data.output.Rdata")###AMD转录组验证集,已取log2
  exprSet.retina.count[1:4,1:4]
  ENSEMBL=ids$ENSEMBL
  RNAseq.retina=as.data.frame(t(exprSet.retina.count[ENSEMBL,]))
  colnames(RNAseq.retina)=ids$gene
  RNAseq.retina=log2(RNAseq.retina+1)
  dim(RNAseq.retina)
  
  RNAseq.retina$disease=as.factor(pdata.retina$group)
  dim(RNAseq.retina)
  str(RNAseq.retina)
  
  test=RNAseq.retina
  
  y=pdata.retina$group
  rf.prob <- predict(final_model, test,type = "prob")
  
  
  rf.ROC = roc(response = y,
               predictor = rf.prob$AMD)
  auc=auc(rf.ROC)
  ##为计算P值
  {
    library(verification)
    y1=ifelse(y=="AMD",1,0)
    ROC.test=roc.area(y1 ,rf.prob$AMD)
    P.value=ROC.test$p.value
  }
  pdf("RNA-ROC_GSE135092_Retina.pdf",useDingbats = F)
  title=paste("GSE135092 Retina samples ROC")
  plot(rf.ROC,
       print.auc=TRUE, print.auc.x=0.4, print.auc.y=0.5,
       # 图像上输出AUC值,坐标为（x，y）
       auc.polygon=TRUE, auc.polygon.col="#fff7f7", # 设置ROC曲线下填充色
       max.auc.polygon=FALSE,  # 填充整个图像
       grid=c(0.5, 0.2), grid.col=c("black", "black"),  # 设置间距为0.1，0.2，线条颜色
       print.thres=TRUE, print.thres.cex=0.9, # 图像上输出最佳截断值，字体缩放倍数
       smooth=F, # 绘制不平滑曲线
       main=title, # 添加标题
       col="#FF2E63",  # 曲线颜色
       legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度
  text(0.4, 0.4, labels=paste("P value =", format.pval(P.value)), col="#FF2E63",adj=c(0, .5)) # 在图上添加P值
  dev.off()
  result.coords <- coords(rf.ROC, "best", best.method="closest.topleft", ret=c("all"))
  result.coords
  # plot(rf.ROC,type = "S",col = "red")
  # text(0.8,0.2, labels = paste0("AUC = ",round(auc,3)))#> [1] 1
  
}