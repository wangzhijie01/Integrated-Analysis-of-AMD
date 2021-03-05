rm(list = ls())   
options(stringsAsFactors = F)


library(RCircos)
load("step7-DEG.DMG.Rdata")
load("gene_location.Rdata")
genelist=rbind(hypo_up,hyper_down)
genelist=merge(genelist,geneid_df,by.x="gene",by.y="gene_name")

genelist = data.frame(Chromosome=paste0("chr",genelist$CHR),
                    chromStart=genelist$start,
                    chromEnd=genelist$end,
                    GeneName=genelist$gene,
                    log2FC=genelist$log2FoldChange,
                    Expr.p=genelist$pvalue,
                    CPG=genelist$cpg,
                    mapinfo=genelist$MAPINFO,
                    deltabeta=genelist$deltaBeta,
                    methy.p=genelist$P.Value
                    )


expression.scatter.data=genelist[!duplicated(genelist$GeneName),]
expression.scatter.data$FC=2^(expression.scatter.data$log2FC)
expression.heatmap.data=expression.scatter.data[,c(1:5)]


gene.label.data=expression.scatter.data[,c(1:4)]
methy.histogram.data=genelist
methy.heatmap.data=data.frame(Chromosome=genelist$Chromosome,
                              chromStart=genelist$mapinfo-50,
                              chromEnd=genelist$mapinfo+50,
                              CPG=genelist$CPG,
                              beta=genelist$deltabeta)
methy.label.data=genelist[,c(1:3,7)]



data(UCSC.HG19.Human.CytoBandIdeogram)
head(UCSC.HG19.Human.CytoBandIdeogram)


# 以下代码使用所有人类染色体表意文字和染色体表意文字内部的10个数据轨道空间来初始化RCircos核心组件。

chr.exclude <- NULL
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,
                            tracks.inside, tracks.outside)

# 核心组件存储在RCircos会话中，每个组件都有一个get方法供高级使用。此外，只需调用函数RCircos.List.Parameters()即可列出所有当前绘图参数。


rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.cyto <- RCircos.Get.Plot.Ideogram()
rcircos.position <- RCircos.Get.Plot.Positions()
RCircos.List.Plot.Parameters()

{
    out.file <- "RCircosDemoHumanGenome.pdf"
    pdf(file=out.file, height=8, width=8,useDingbats = F)
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()
    
    name.col <- 4
    side <- "in"
    track.num <- 1
    RCircos.Gene.Connector.Plot(gene.label.data,
                                track.num, side)###基因名横条
    track.num <- 2
    RCircos.Gene.Name.Plot(gene.label.data,
                           name.col,track.num, side)###基因名label
    
    
    data.col <- 5
    track.num <- 4
    side <- "in"
    RCircos.Heatmap.Plot(expression.heatmap.data, data.col,min.value=-2,max.value = 2,
                         track.num, side)
    
    data.col <- 5;
    track.num <- 5;
    side <- "in";
    RCircos.Heatmap.Plot(methy.heatmap.data, data.col, track.num, ,min.value=-0.1,max.value = 0.1,side);
    dev.off()
}


##pie 
slices = c(15,24,11,19)
lbls = c("hyper_down","hyper_up","hypo_down","hypo_up")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
pdf("pie.pdf",useDingbats = F)
pie(slices, labels = lbls, 
    main="Pie Chart of Species\n (with sample sizes)")
dev.off()
    