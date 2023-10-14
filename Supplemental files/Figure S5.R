library(stats)
library(dplyr)
library(tidyr)
load("Cel_dup_N4Cel.RData")
load("celltype.RData")
load("sc_exp.RData")

for(i in 1:ncol(sc_exp)){
  cur = filter(celltype, index == colnames(sc_exp)[i])
  colnames(sc_exp)[i] = cur$celltype2
}

for(i in 1:nrow(sc_exp)){
  sc_exp[i,"max_cell"] = colnames(sc_exp)[which(sc_exp[i,] == max(as.numeric(sc_exp[i,]),na.rm = T))[[1]]]
  if(i %% 5000 ==0){
    cat(i,"rows finished \n")
  }
}
n4cel_sc = intersect(unique(Cel_dup_N4Cel$WBGene),unique(rownames(sc_exp)))
n4cel_sc_others = setdiff(n4cel_sc,rownames(n4cel_cls10))
sc_exp_others = sc_exp[n4cel_sc_others,]

cell1_16 = celltype[which(celltype[,"celltype2"] == "P0_1-cell") : nrow(celltype),]
tpm_lin = select(sc_exp,intersect(celltype$celltype2,c(cell1_16$celltype2,lineage_allcells$X)))
tpm_emb = select(sc_exp,intersect(celltype$celltype2,c(emb_late$type)))
tpm_l2 = select(sc_exp,intersect(celltype$celltype2,c(L2_cell$celltype)))
tpm_l4 = select(sc_exp,intersect(celltype$celltype2,c(L4_cell$celltype)))

tpm_all = list(tpm_lin,tpm_emb,tpm_l2,tpm_l4)
type4 = c("lineage","late_emb","l2","l4")
tpm_all_others_exp = list()
tpm_all_others_z = list()
for(i in 1:length(type4)){
  cur = log2(tpm_all[[i]] +1)
  cur2 = cur[rownames(sc_exp_others),]
  cur2_mean = apply(cur2,2,mean2)
  cur2_mean = data.frame(cur2_mean)
  colnames(cur2_mean) = "mean_exp"
  cur2_mean$type = type4[i]
  tpm_all_others_exp[[i]] = cur2_mean
  
  cur2_z = t(apply(cur2,1,zscore))
  cur2_z_mean = apply(cur2_z,2,mean2)
  cur2_z_mean = data.frame(cur2_z_mean)
  colnames(cur2_z_mean) = "mean_zscore"
  cur2_z_mean$type = type4[i]
  tpm_all_others_z[[i]] = cur2_z_mean
}
tpm_all_others_exp2 = do.call(rbind,tpm_all_others_exp)
tpm_all_others_z2 = do.call(rbind,tpm_all_others_z)
tpm_all_others_exp2$celltype = rownames(tpm_all_others_exp2)
tpm_all_others_z2$celltype = rownames(tpm_all_others_z2)
tpm_all_others_exp2_zscore2 = merge(tpm_all_others_exp2,tpm_all_others_z2,by = "celltype")


