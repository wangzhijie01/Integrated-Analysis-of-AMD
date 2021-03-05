
rm(list = ls())   
options(stringsAsFactors = F)

load("step3.1-sub_deg.Rdata")
load("D:/bioinformation/R/甲基化/AMD/转录组/anno_nrDEG_macular.Rdata")


{
  
  DMG.DEG=merge(sub_deg,nrDEG,by.x="gene",by.y="SYMBOL")
  table(DMG.DEG$change.x,DMG.DEG$change.y)
  hypo_up=DMG.DEG[DMG.DEG$change.x=="DOWN"&DMG.DEG$change.y=="UP",]
  hypo_DOWN=DMG.DEG[DMG.DEG$change.x=="DOWN"&DMG.DEG$change.y=="DOWN",]
  
  hyper_down=DMG.DEG[DMG.DEG$change.x=="UP"&DMG.DEG$change.y=="DOWN",]
  hyper_up=DMG.DEG[DMG.DEG$change.x=="UP"&DMG.DEG$change.y=="UP",]
  
}
write.csv(hyper_down,file = "hyper_down.csv")
write.csv(hypo_up,file = "hypo_up.csv")
save.image("step7-DEG.DMG.Rdata")



# 画维恩图
{
  library(venn)
  library(RColorBrewer)

  #读入作图文件
  sub_deg=sub_deg[!sub_deg$feature=="IGR",]
  feature=unique(sub_deg$feature)
  hypo=as.character(unique(sub_deg$gene[sub_deg$change=="DOWN"]))
  hyper=as.character(unique(sub_deg$gene[sub_deg$change=="UP"]))
  gene_down=nrDEG$SYMBOL[nrDEG$change=="DOWN"]
  gene_up=nrDEG$SYMBOL[nrDEG$change=="UP"]
  venn_list1=list(Hypomethylated=hypo,Gene_Down=gene_down,Gene_Up=gene_up)
  venn_list2=list(Hypermethylated=hyper,Gene_Down=gene_down,Gene_Up=gene_up)
  display.brewer.pal(3, "Set2")
  mycolors<-brewer.pal(3, "Set2")
  #作图，详情使用 ?venn 查看帮助（注：它的作图函数也是 venn()，和上述 gplots 包存在冲突）
  pdf("Hypo-DMGs&DEG.pdf",width = 6,height = 6,useDingbats = F)
  venn(venn_list1, opacity=0.5,borders=T,box=F,ilcs = 0.8,zcolor = mycolors)
  dev.off()
  
  pdf("Hyper-DMGs&DEG.pdf",width = 6,height = 6,useDingbats = F)
  venn(venn_list2, opacity=0.5,borders=T,box=F,ilcs = 0.8,zcolor = mycolors)
  dev.off()
  
}
