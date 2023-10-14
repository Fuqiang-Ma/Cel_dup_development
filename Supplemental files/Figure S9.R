library(stats)
library(purrr)
library(dplyr)
load("bulk_dev_RNAseq.RData")
load("stress_RNAseq.RData")
load("Cel_dup_N4Cel.RData")
load("cluster_1_2_3.RData")
load("sc_exp.RData")
load("celltype.RData")
load("L4_exp.RData")

zscore = function(x){
  x = (x - mean(x))/sd(x)
}

mean2 = function(x){
  mean(x,na.rm = T)
}

bulk_RNAseq_n4cel = bulk_RNAseq[bulk_RNAseq$WBGene %in% unique(Cel_dup_N4Cel$WBGene),]
n4cel_null = filter(bulk_RNAseq_n4cel, max <= 1e-10)
n4cel_low = filter(bulk_RNAseq_n4cel, max <= 10 & max > 1e-10)
n4cel_medium = filter(bulk_RNAseq_n4cel, max > 10 & max <= 100)
n4cel_high = filter(bulk_RNAseq_n4cel, max > 100)

bulk_RNAseq_no_n4cel = bulk_RNAseq[bulk_RNAseq$WBGene %in% unique(setdiff(bulk_RNAseq$WBGene,Cel_dup_N4Cel$WBGene)),]
genes_null2 = filter(bulk_RNAseq_no_n4cel, max <= 1e-10)
genes_low2 = filter(bulk_RNAseq_no_n4cel, max <= 10 & max > 1e-10)
genes_medium2 = filter(bulk_RNAseq_no_n4cel, max > 10 & max <= 100)
genes_high2 = filter(bulk_RNAseq_no_n4cel, max > 100)

rec_z_cls2$WBGene = rownames(rec_z_cls2)
cls1 = filter(rec_z_cls2, cls == 1)
cls2 = filter(rec_z_cls2, cls == 2)
cls3 = filter(rec_z_cls2, cls == 3)

stress = filter(stress, condition == c("heatshock","Bortezomib","Tunicamycin","cold1","hypoxia","osmotic2","starvation2"))
###genes enrichment with fc > 3
stress3 = filter(stress,stress[,"log2FoldChange"] > 1.58)
Cel_dup_N4Cel$type = "recent_dup"
genes_1to1$type = "single_copy"
type = c("n4cel_null","n4cel_low","n4cel_medium","n4cel_high","cls1","cls2","cls3","other_genes_0","other_genes_0_10","other_genes_10_100","other_genes_100","single_copy_genes")
exp_list = list(n4cel_null,n4cel_low,n4cel_medium,n4cel_high,cls1,cls2,cls3,genes_null2,genes_low2,genes_medium2,genes_high2,genes_1to1)
exp_enr = data.frame(type)
library(stats)
for(i in 1:length(type)){
  cur = exp_list[[i]][exp_list[[i]]$WBGene %in% unique(stress3$WBGene),]
  exp_enr[i,"FC2"] = nrow(cur)
  exp_enr[i,"total"] = nrow(exp_list[[i]])
  exp_enr[i,"enrich"] = phyper(c(exp_enr[i,"FC2"]-1),exp_enr[i,"total"],c(19997 - exp_enr[i,"total"]),length(unique(stress3$WBGene)),lower.tail = F)
  exp_enr[i,"enrich2"] = -log10(exp_enr[i,"enrich"])
}

###plasticity
stress3 = filter(stress,log2FoldChange>0)
type = c("n4cel","n4cel_null","n4cel_low","n4cel_medium","n4cel_high","cls1","cls2","cls3","other_genes_0","other_genes_0_10","other_genes_10_100","other_genes_100","single_copy_genes")
exp_list = list(Cel_dup_N4Cel,n4cel_null,n4cel_low,n4cel_medium,n4cel_high,cls1,cls2,cls3,genes_null2,genes_low2,genes_medium2,genes_high2,genes_1to1)
exp_list2 = list()
for(i in 1:length(type)){
  cur = stress3[stress3$WBGene %in% unique(exp_list[[i]]$WBGene),]
  colnames(cur) = c("WBGene",type[i])
  exp_list2[[i]] = cur
}
exp_list2_3 = Reduce(function(x,y) merge(x,y,by = "WBGene",all=T),exp_list2)

###versatility
stress3 = filter(stress,stress[,"log2FoldChange"] > 1.58)
type = c("n4cel","n4cel_null","n4cel_low","n4cel_medium","n4cel_high","cls1","cls2","cls3","other_genes_0","other_genes_0_10","other_genes_10_100","other_genes_100","single_copy_genes")
exp_list = list(Cel_dup_N4Cel,n4cel_null,n4cel_low,n4cel_medium,n4cel_high,cls1,cls2,cls3,genes_null2,genes_low2,genes_medium2,genes_high2,genes_1to1)
exp_list_versat = list()
exp_list_versat_mean = c()
for(i in 1:length(type)){
  cur = stress3[stress3$WBGene %in% unique(exp_list[[i]]$WBGene),]
  cur_sum = data.frame(table(cur$WBGene))
  colnames(cur_sum) = c("WBGene","Freq")
  cur_sum = filter(cur_sum,Freq >0)
  cur_sum$type = type[i]
  exp_list_versat[[i]] = cur_sum
  exp_list_versat_mean = append(exp_list_versat_mean,mean(cur_sum$Freq))
}
exp_list_versat_2 = do.call(rbind,exp_list_versat)
exp_list_versat_mean2 = data.frame(type,exp_list_versat_mean)

### L4 enrichment_cls3
aaa2 = celltype[celltype$celltype2 %in% colnames(L4_exp),]
L4_tpm = sc_exp[,aaa2$index]
for(i in 1:ncol(L4_tpm)){
  for(j in 1:nrow(aaa2)){
    if(colnames(L4_tpm)[i] == aaa2[j,"index"]){
      colnames(L4_tpm)[i] = aaa2[j,"celltype2"]
    }
  }}
L4_zscore = data.frame(t(apply(log2(L4_tpm + 1),1,zscore)))
upregulate = list()
for(i in 1:length(type4)){
  cur = stress3[stress3$WBGene %in% unique(exp_list2[[i]]$WBGene),]
  cur_L4 = L4_zscore[rownames(L4_zscore) %in% unique(cur$WBGene),]
  cur2 = filter(cur_L4,!is.na(cur_L4[,"AWB"]))
  cur2_mean = data.frame(sort(apply(cur2,2,mean2),decreasing = T))
  colnames(cur2_mean) = "L4"
  cur2_mean$type = type4[i]
  cur2_mean$celltype = rownames(cur2_mean)
  upregulate[[i]] = cur2_mean
}
upregulate2 = do.call(rbind,upregulate)


