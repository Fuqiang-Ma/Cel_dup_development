library(dplyr)
load("cell_emb_late.RData")
load("cell_L4.RData")
load("sc_exp.RData")
load("celltype.RData")
load("cluster_1_2_3.RData")

zscore = function(x){
  x = (x - mean(x))/sd(x)
}

###all cells at 440 mins enrichment
min440 = filter(emb_late,Start_time < 440 & End_time >= 440)
names(min440)[names(min440) == "type"] = "Precursor2"
cell1 = sort(unique(min440$Precursor2))
cell2 = celltype2_3[celltype2_3$celltype2 %in% cell1,]
celltype_440_tpm = TPM_mean2[,cell2$index]
for(i in 1:ncol(celltype_440_tpm)){
  for(j in 1:nrow(cell2)){
    if(colnames(celltype_440_tpm)[i] == cell2[j,"index"]){
      colnames(celltype_440_tpm)[i] = cell2[j,"celltype2"]
    }
}}

celltype_440_zscore = data.frame(t(apply(log2(celltype_440_tpm + 1),1,zscore)))
celltype_440_exp = log2(celltype_440_tpm + 1)
cell440 = list(celltype_440_zscore,celltype_440_exp)
type = c("min440_zscore","min440_exp")
cls2 = filter(rec_z_cls2,cls == 2)
for(i in 1:length(cell440)){
  cur = cell440[[i]][rownames(cls2),]
  cur = filter(cur,!is.na(cur[,5]))
  cur_mean = data.frame(sort(apply(cur,2,mean2)))
  colnames(cur_mean) = "min455"
  write.csv(cur_mean,file = paste0(type[i],".csv"),row.names = T)
}

### cells at L4 enrichment
cell1 = sort(unique(L4_cell$celltype))
cell2 = celltype2_3[celltype2_3$celltype2 %in% cell1,]
L4_tpm = TPM_mean2[,cell2$index]
for(i in 1:ncol(L4_tpm)){
  for(j in 1:nrow(cell2)){
    if(colnames(L4_tpm)[i] == cell2[j,"index"]){
      colnames(L4_tpm)[i] = cell2[j,"celltype2"]
    }
}}

L4_zscore = data.frame(t(apply(log2(L4_tpm + 1),1,zscore)))
L4_exp = log2(L4_tpm + 1)
L4 = list(L4_zscore,L4_exp)
type = c("L4_zscore","L4_exp")
cls3 = filter(rec_z_cls2,cls == 3)
for(i in 1:length(L4)){
  cur = L4[[i]][rownames(cls3),]
  cur = filter(cur,!is.na(cur[,5]))
  cur_mean = data.frame(sort(apply(cur,2,mean2)))
  colnames(cur_mean) = "L4"
  write.csv(cur_mean,file = paste0(type[i],".csv"),row.names = T)
}

