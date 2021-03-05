rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load("anno_nrDEG_macular.Rdata")

GO_diff=read.csv("macular_gene_up_GO_all.csv")
GO_diff=GO_diff[,c(2:4,10,8)]
colnames(GO_diff)=c("Category","ID","Term","Genes","adj_pvalue")
GO_diff$Genes=gsub("/",", ",GO_diff$Genes)

KEGG_up=read.csv("kegg_macular__kk.up.csv")
KEGG_up=data.frame(Category="KEGG",KEGG_up[,c(2,3,9,7)])
colnames(KEGG_up)=c("Category","ID","Term","Genes","adj_pvalue")
KEGG_up$Genes=gsub("/",", ",KEGG_up$Genes)

GO_KEGG=rbind(GO_diff,KEGG_up)

genelist=nrDEG[,c(5,2,3,6)]
colnames(genelist)=c("ID","logFC","adj.P.Val","ENTREZID")
##需要呈现的terms
process=c("IL-17 signaling pathway","TNF signaling pathway",
          "ECM-receptor interaction","Cytokine-cytokine receptor interaction","extracellular structure organization")

genes=genelist[,c(1:2)]
colnames(genes)=c("ID","logFC")
# row=sample(1:439,100)
# genes=genes[row,]



library(GOplot)


EC=list()
EC$david=GO_KEGG
EC$genelist=genelist
EC$process=process
EC$genes=genes


# Generate the plotting object
circ <- circle_dat(EC$david, EC$genelist)

# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))

# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')

#Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)

# Generate a circular visualization of selected terms
IDs <- c("GO:0043062" ,"GO:0030198", "GO:0051302", "GO:0045907")
GOCircle(circ, nsub = IDs)


# Define a list of genes which you think are interesting to look at. The item EC$genes of the toy 
# sample contains the data frame of selected genes and their logFC. Have a look...
tmp=(EC$genes)

EC$process
# Now it is time to generate the binary matrix
chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)


# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = EC$genes)
# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = EC$process)

# Create the plot
chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)

pdf("chordplot.pdf",width = 12.5,height = 12.5,useDingbats = F)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
dev.off()


pdf("chordplot_kegg_UP_5.pdf",width = 12.5,height = 12.5,useDingbats = F)
KEGG_terms=c("IL-17 signaling pathway","Neuroactive ligand-receptor interaction","TNF signaling pathway",
             "ECM-receptor interaction","NF-kappa B signaling pathway")
chord <- chord_dat(data = circ, genes = EC$genes, process = KEGG_terms)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
dev.off()


pdf("chordplot_GO_UP_5.pdf",width = 12.5,height = 12.5,useDingbats = F)
GO_terms=c("extracellular structure organization","extracellular matrix organization","neutrophil chemotaxis",
           "neutrophil migration","cytokine activity","chemokine activity")
chord <- chord_dat(data = circ, genes = EC$genes, process = GO_terms)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
dev.off()


# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(2, 0), gene.order = 'logFC')

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-8], nlfc = 0)


# Now we create the heatmap with logFC values and user-defined colour scale
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))


GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)

GOCluster(circ, EC$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))

