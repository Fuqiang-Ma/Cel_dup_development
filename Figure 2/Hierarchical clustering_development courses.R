library(dplyr)
library(stats)
load("bulk_dev_RNAseq.RData")
load("cluster_1_2_3.RData")

zscore = function(x){
  x = (x - mean(x))/sd(x)
}

for(i in 1:nrow(bulk_RNAseq)){
  bulk_RNAseq[i,"max"] = max(as.numeric(bulk_RNAseq[i,2:49]))
}

bulk_RNAseq2 = filter(bulk_RNAseq,max>10)
rownames(bulk_RNAseq2) = bulk_RNAseq2$WBGene
bulk_RNAseq2 = select(bulk_RNAseq2,-c("WBGene","max"))
bulk_RNAseq_z = t(apply(log2(bulk_RNAseq2+1),1,zscore))
#clustering
rec_z = bulk_RNAseq_z[rownames(bulk_RNAseq_z) %in% unique(Cel_dup_N4Cel$WBGene),]
rec_z[is.na(rec_z)] = 0
rec_z_cls = hclust(dist(rec_z), method = "complete")
rec_z_cls2 = data.frame(cutree(rec_z_cls,3))
colnames(rec_z_cls2) = "cls"

for(i in 1:nrow(rec_z_cls2)){
   if(rec_z_cls2[i,"cls"] == 2){
     rec_z_cls2[i,"cls"] =3
     } else if(rec_z_cls2[i,"cls"] == 3){
       rec_z_cls2[i,"cls"] =2
   }
}

#zscore
cls1_2_3_z = data.frame(matrix(nrow = 3, ncol = ncol(bulk_RNAseq_z)))
rownames(cls1_2_3_z) = paste0("cluster",1:3)
colnames(cls1_2_3_z) = colnames(bulk_RNAseq_z)
for(j in 1:nrow(cls1_2_3_z)){
  cur = filter(rec_z_cls2,cls == j)
  cur2 = bulk_RNAseq_z[rownames(bulk_RNAseq_z) %in% rownames(cur),]
  cls1_2_3_z[j,] = apply(cur2,2,mean)
}

write.csv(cls1_2_3_z,file = "n4cel_cls10_exp.csv",row.names=T)

#expression
bulk_RNAseq2 = log2(bulk_RNAseq2+1)
cls1_2_3_exp = data.frame(matrix(nrow = 3, ncol = ncol(bulk_RNAseq2)))
rownames(cls1_2_3_exp) = paste0("cluster",1:3)
colnames(cls1_2_3_exp) = colnames(bulk_RNAseq2)
for(j in 1:nrow(cls1_2_3_exp)){
  cur = filter(rec_z_cls2,cls == j)
  cur2 = bulk_RNAseq2[rownames(bulk_RNAseq2) %in% rownames(cur),]
  cls1_2_3_exp[j,] = apply(cur2,2,mean)
}

write.csv(cls1_2_3_exp,file = "n4cel_cls10_exp.csv",row.names=T)


