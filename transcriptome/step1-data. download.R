rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
# 注意查看下载文件的大小，检查数据 
f='GSE135092_eSet.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE135092', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE135092_eSet.Rdata')  ## 载入数据
class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])
gset
# assayData: 44149 features, 41 samples

# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]] #

# GPL13667
dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
boxplot(dat[,1:10],las=2)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData
## 挑选一些感兴趣的临床表型。
pdata.raw=pd[,c(2,33,35)]
colnames(pdata.raw)=c("geo_accession","group","tissue")
head(pdata.raw)
#筛选RPE样本
pdata=pdata.raw[grep("RPE",pdata.raw$tissue),]
dim(pdata)
table(pdata$tissue,pdata$group)

pdata.retina=pdata.raw[grep("Retina",pdata.raw$tissue),]


##read raw count
# GPL570
if (F) {
  files = dir("GSE135092_RAW/")
  ##查看多少个文件
  length(files)
  files=files[substr(files,1,10)%in%pdata$geo_accession]
  length(files)
  ### 551
  ## View(files)可以看一下files具体啥样
  ##先读取rawdata文件夹下第一个文件夹中的数据看对不对，如果对的话就采用for循环读取所有数据，并将所有数据合并在一起
  ##paste0("rawdata","/",files[1],"/",dir(paste0("rawdata","/",files[1]))[1] )即为"rawdata/0100d06e-cc9d-4ddc-a9b2-93bfea759f7a/64bfb99a-508d-4f1e-9365-58745a59d518.htseq.counts.gz"，这样是不是就比较容易理解呀！
  out = read.table(gzfile(paste0("GSE135092_RAW/",files[1])),header = T)
  ##把读取的数据的第二列即表达值的那列改为文件名（文件名代表样本名）
  colnames(out)[2:3] <- c(paste0(substr(files[1],1,10),"-RPKM"),paste0(substr(files[1],1,10),"-count"))
  ##采用for循环读取所有数据，并将所有数据合并在一起
  for(x in c(2:length(files))){
    temp = read.table(gzfile(paste0("GSE135092_RAW/",files[x])),header = T)##读取数据
    colnames(temp)[2:3] <- c(paste0(substr(files[x],1,10),"-RPKM"),paste0(substr(files[x],1,10),"-count"))##把读取的数据的第二列即表达值的那列改为文件名
    out <- merge(out,temp,by="ID_REF")##合并所有数据
    print(paste0("GSE135092_RAW/",files[x]))
  }

  rownames(out)=out[,1]
  out=out[,-1]
  exprSet.RPKM=out[,grep("RPKM",colnames(out))]
  exprSet.count=out[,grep("count",colnames(out))]
  colnames(exprSet.count)=substr( colnames(exprSet.count),1,10)
  colnames(exprSet.RPKM)=substr( colnames(exprSet.RPKM),1,10)
}

save(exprSet.RPKM,exprSet.count,pdata,file = 'step1-output.Rdata')

##read raw count
# GPL570
if (F) {
  files = dir("GSE135092_RAW/")
  ##查看多少个文件
  length(files)
  files=files[substr(files,1,10)%in%pdata.retina$geo_accession]
  length(files)
  ### 551
  ## View(files)可以看一下files具体啥样
  ##先读取rawdata文件夹下第一个文件夹中的数据看对不对，如果对的话就采用for循环读取所有数据，并将所有数据合并在一起
  ##paste0("rawdata","/",files[1],"/",dir(paste0("rawdata","/",files[1]))[1] )即为"rawdata/0100d06e-cc9d-4ddc-a9b2-93bfea759f7a/64bfb99a-508d-4f1e-9365-58745a59d518.htseq.counts.gz"，这样是不是就比较容易理解呀！
  out = read.table(gzfile(paste0("GSE135092_RAW/",files[1])),header = T)
  ##把读取的数据的第二列即表达值的那列改为文件名（文件名代表样本名）
  colnames(out)[2:3] <- c(paste0(substr(files[1],1,10),"-RPKM"),paste0(substr(files[1],1,10),"-count"))
  ##采用for循环读取所有数据，并将所有数据合并在一起
  for(x in c(2:length(files))){
    temp = read.table(gzfile(paste0("GSE135092_RAW/",files[x])),header = T)##读取数据
    colnames(temp)[2:3] <- c(paste0(substr(files[x],1,10),"-RPKM"),paste0(substr(files[x],1,10),"-count"))##把读取的数据的第二列即表达值的那列改为文件名
    out <- merge(out,temp,by="ID_REF")##合并所有数据
    print(paste0("GSE135092_RAW/",files[x]))
  }
  
  rownames(out)=out[,1]
  out=out[,-1]
  exprSet.retina.RPKM=out[,grep("RPKM",colnames(out))]
  exprSet.retina.count=out[,grep("count",colnames(out))]
  colnames(exprSet.retina.count)=substr( colnames(exprSet.retina.count),1,10)
  colnames(exprSet.retina.RPKM)=substr( colnames(exprSet.retina.RPKM),1,10)
}

save(exprSet.retina.RPKM,exprSet.retina.count,pdata.retina,file = 'step1-retina.data.output.Rdata')

