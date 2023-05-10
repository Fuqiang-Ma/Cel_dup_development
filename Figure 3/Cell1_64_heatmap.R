library(dplyr)
library(RColorBrewer)
library(gplots)
load("sc_exp.RData")
load("celltype.RData")
load("cell1_64.RData")
load("cell1_64_2.RData")
load("Single_copy_genes.RData")
load("cluster_1_2_3.RData")

mean2 = function(x){
  mean(x,na.rm = T)
}

cell_list = list()
for(i in 1:nrow(cell1_64)){
  cur = cell1_64_2[cell1_64_2$celltype3 %in%  as.character(cell1_64[i,]),]
  cell = data.frame(as.character(cell1_64[i,]),1:length(as.character(cell1_64[i,])))
  colnames(cell) = c("cel","order")
  cur2 = merge(cell, cur,by.x = "cel",by.y = "celltype3",all.x = T)
  cur2 = arrange(cur2,cur2$order)
  
  cell2 = data.frame(matrix(nrow = nrow(TPM_mean2),ncol = nrow(cur2)))
  rownames(cell2) = rownames(TPM_mean2)
  colnames(cell2) = cur2$cel
  #cell2 = TPM_mean2[,cur2$index]
  for(j in 1:ncol(cell2)){
    cur3 = filter(cur2,cel == colnames(cell2)[j])
    if(!is.na(cur3[1,"index"])){
      cell2[,j] = TPM_mean2[,cur3$index]
    }
  }
  cell_list[[i]] = log2(cell2+1) #note exp value has been log2 transformed
}

### heatmap
cls1 = filter(rec_z_cls2,cls ==1)
rownames(genes_1to1) = genes_1to1$WBGene
genes = list(cls1,genes_1to1)
name = c("cls1","single_copy")

for(i in 1:length(genes)){
  cell_rec = data.frame(matrix(nrow = nrow(cell1_64), ncol = 7))
  cell1_64$V9 = cell1_64$V8
  cell1_64[c(7,17,21,22,25,26,28,29),"V9"] = c("Exxx_Ep","ABalppp/ABpraaa_ABa","ABarpap/ABplaaa_ABa","ABarpax","ABpraxx","ABprpax","ABplaxx","ABplpax")
  rownames(cell_rec) = cell1_64$V9
  colnames(cell_rec) = c("C1","C2","C4","C8","C16","C32","C64")

for(j in 1:nrow(cell_rec)){
  cur = cell_list[[j]][rownames(cell_list[[j]]) %in% rownames(genes[[i]]),]
  cell_rec[j,] = apply(cur,2,mean2) # use mean or median accordingly
}

my_palette <- colorRampPalette(c("lightgreen", "black", "red"))(n = 299)
col_breaks = c(seq(2,3.331,length=100),  # for green
               seq(3.3311,4.631,length=100),           # for black
               seq(4.6311,6,length=100))
pdf(paste0("emb1_64_",name[i],".pdf"),height = 15,width=25)
heatmap.2(as.matrix(cell_rec),
          density.info="none",  # turns off density plot inside color legend
          key=F, keysize=1.0,  ## change F as T if want color scaling key
          scale="none",
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits # use this if want to limit color range
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
          srtCol=22,            # turn off column clustering
          margins =c(11,11),
)
dev.off()
}

