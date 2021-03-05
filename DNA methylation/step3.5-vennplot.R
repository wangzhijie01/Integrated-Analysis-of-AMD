rm(list = ls())
options(stringsAsFactors = F)

load(file="step3.1-sub_deg.Rdata")
library(venn)
library(RColorBrewer)

#读入作图文件
sub_deg=sub_deg[!sub_deg$feature=="IGR",]
feature=unique(sub_deg$feature)
venn_list=list()
for (i in feature) {
  sub=sub_deg[sub_deg$feature==i,]
  venn_list[[i]]=as.character(sub$gene)
}
display.brewer.pal(6, "Set3")
mycolors<-brewer.pal(6, "Set3")
#作图，详情使用 ?venn 查看帮助（注：它的作图函数也是 venn()，和上述 gplots 包存在冲突）
pdf("DMGs_Venn.pdf",width = 6,height = 6,useDingbats = F)
venn(venn_list, opacity=0.5,borders=T,ilcs = 0.8,zcolor = mycolors)
dev.off()
