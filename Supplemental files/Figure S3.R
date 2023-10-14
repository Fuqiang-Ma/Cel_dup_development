library(dplyr)
library(RColorBrewer)
library(gplots)
load("cluster_1_2_3.RData")
load("bulk_dev_RNAseq.RData")

rownames(bulk_RNAseq) = bulk_RNAseq$WBGene
bulk_RNAseq = select(bulk_RNAseq,-c("WBGene","max"))
my_palette <- colorRampPalette(c("lightgreen", "black", "red"))(n = 299)
col_breaks = c(seq(-6,-2.01,length=100),  # for green
               seq(-2,2,length=100),           # for black
               seq(2.01,6,length=100))
for(j in 1:3){
  cur = filter(rec_z_cls2,cls == j)
  cur2 = bulk_RNAseq[rownames(bulk_RNAseq) %in% rownames(cur),]
  cur3 = as.matrix(cur2)
  pdf(paste0("cls_",j,"_zscore.pdf",sep=""),height = 15,width=25)
  heatmap.2(cur3,
          density.info="none",
          key=F, keysize=1.0,
          scale="none",
          trace="none",            
          col=my_palette,       
          breaks=col_breaks,    
          dendrogram="row",
          Colv="NA",
          labRow = "none",
          srtCol=22,
          margins =c(11,11))
dev.off()
}


