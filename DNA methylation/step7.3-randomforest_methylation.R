
rm(list = ls())   
options(stringsAsFactors = F)

load("step7-DEG.DMG.Rdata")##merge
load("step2-champ_myNorm.Rdata")###甲基化矩阵
methylation.pdata=myLoad$pd

###数据预处理###

if (F) {
  CPG=c(hypo_up$cpg,hyper_down$cpg)
  CPG=unique(CPG)
  gene=unique(c(hypo_up$gene,hyper_down$gene))
  
  methylation=as.data.frame(t(myNorm2[CPG,]))
  str(methylation)
  methylation$disease=as.factor(methylation.pdata$Sample_Group)
  dim(methylation)
  str(methylation)
}


##randome-forest 建模##

library(randomForest)

##参数筛选##
if(T){
  # 筛选ntree
  set.seed(100)
  rf_ntree <- randomForest(disease~.,data=methylation,ntree=1000)
  pdf("methylation_ntree_chosen.pdf",useDingbats = F)
  plot(rf_ntree)
  dev.off()
  # 筛选mtry
  n <- ncol(methylation)-1
  errRate <- c()
  
  # 从图中可以看到，当ntree=400时，模型内的误差就基本稳定了，出于更保险的考虑，我们确定ntree值为100。
  set.seed(100)
  for (i in 1:n){ 
    
    m <- randomForest(disease~.,data=methylation,mtry=i,proximity=TRUE) 
    
    err<-mean(m$err.rate[,1])
    
    errRate[i] <- err
  }  
  print(errRate)
  m= which.min(errRate)  
  print(m)
  pdf("methylation_ntry_chosen.pdf",useDingbats = F)
  matplot(1:n , errRate, pch=19 , col=c("red"),type="p",ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
  legend("topright",legend=c("Out of Bag Error"),pch=19, col=c("red"))
  dev.off()
  #选择最优mtry参数值
  
}

##importance 筛选基因排序##
if (F) {
  set.seed(100)
  rf_output=randomForest(disease~., data =methylation ,importance = TRUE,mtry=2,ntree=200,
                         proximity=TRUE)
  rf_importances=importance(rf_output, scale=FALSE)
  choose_gene=rownames(tail(rf_importances[,],10))
  head(rf_importances)
  gene.order=rownames(rf_importances)[order(rf_importances[,4],decreasing = T)]
  write.csv(rf_importances,file = "methylation_importance.csv")
  
  pdf("methylation-rf_importance.pdf",useDingbats = F)
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
  grid=expand.grid(mtry=c(2))
  # train the model
  model_RNA <- train(disease~., data =methylation , trControl=train_control,method="rf",tuneGrid =grid,ntree=200)
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
    tmp.data = methylation[,c(gen.enroll,"disease")]
    #Leave One Out Cross Validation
    train_control <- trainControl(method="LOOCV",
                                  summaryFunction=twoClassSummary, classProbs=TRUE,##把评价改成ROC，不然是accuracy
                                  savePredictions = T)
    
    #Repeated k-fold Cross Validation
    # train_control <- trainControl(method="repeatedcv", number=10, repeats=3,
    #                               summaryFunction=twoClassSummary, classProbs=TRUE,##把评价改成ROC，不然是accuracy
    #                               savePredictions = T)##把predict值存储下来画ROC用
    grid = expand.grid(mtry=2)
    # train the model
    model_RNA <- train(disease~., data =tmp.data,trControl=train_control,method="rf",tuneGrid =grid,ntree=200)
    # summarize results
    print(model_RNA)
    ROC[i] = model_RNA$results$ROC[1]
    RNA_model[[i]] = model_RNA
  }
  
  
  ROC
  
  pdf("Gene methylation classifier_chosen_gene.pdf",useDingbats = F)
  plot(1:n,ROC, pch=19 , col=c("black"),type="p",cex=1,
       ylab="AUC of ROC curve",xlab="Number of Predictors in classifier",
       main="Gene methylation classifier based on CpGs")
  dev.off()
  
  ntmp=which.max(ROC)
  ntmp
  final_model = RNA_model[[ntmp]]
  ROC[ntmp]
  save(ROC,RNA_model,final_model,file = "step7.3-methylation-final_rf_model.Rdata")
  
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
  pdf("methylation_ROC_train.pdf",useDingbats = F)
  title=paste("Top",ntmp,"CpGs methylation classifier")
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
}



