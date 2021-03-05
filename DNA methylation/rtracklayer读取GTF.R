library(rtracklayer)
gtf1 <- rtracklayer::import('D:/bioinformation/GTF注释文件/Homo_sapiens.GRCh38.94.gtf')
gtf_df <- as.data.frame(gtf1)
test <- gtf_df[1:20,]
View(test)


geneid_df <- dplyr::select(gtf_df,c(gene_name,seqnames,start,end,type))
geneid_df = geneid_df[geneid_df$type=="gene",]

save(geneid_df,file="gene_location.Rdata")
