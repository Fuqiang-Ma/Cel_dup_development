library(stats)
library(purrr)
library(dplyr)
load("bulk_dev_RNAseq.RData")
load("stress_RNAseq.RData")
load("Cel_dup_N4Cel.RData")
load("Single_copy_genes.RData")
load("cluster_1_2_3.RData")
load("domain_elegans.RData")
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

###RNAseq FC3 volcano plot
stress = filter(stress, condition %in% c("Db10","BT","PA14","Nparisii_DAVID","orsay_DAVID","Pathogen_Enterococcus_faecalis","Pathogen_Penicillium_brevicompactum","Myzocytiopsis_humicola","Pseudomonas_vranovensis"))
stress3 = filter(stress,abs(stress[,"log2FoldChange"]) > 1.58)
NHR = domain_elegans[grep("\\bHormone_recep|zf-C4\\b",domain_elegans$domain),]
FBP = domain_elegans[grep("\\bF-box|\\bFTH\\b|\\bHTH_48\\b|\\bFBA_2\\b",domain_elegans$domain),]
MATH = domain_elegans[grep("\\bMATH\\b",domain_elegans$domain),]
LEC = domain_elegans[grep("\\bLectin_C\\b",domain_elegans$domain),]
domain_list = list(NHR,FBP,MATH,LEC)
domain_list2 = list()
type3 = c("NHR","FBP","MATH","LEC") #this is the same order with domain_list
for(i in 1:length(type3)){
  cur = stress3[stress3$WBGene %in% unique(domain_list[[i]]$Wb),]
  cur$type = type3[i]
  domain_list2[[i]] = cur
}
domain_list2_2= do.call(rbind,domain_list2)
others = stress3[stress3$WBGene %in% setdiff(stress3$WBGene, domain_list2$WBGene),]
others$type = "others"
domain_list2_3 = rbind(domain_list2_2,others)
domain_list2_3$pvalue2 = -log10(domain_list2_3$pvalue)

###genes enrichment with fc > 3
stress3 = filter(stress,stress[,"log2FoldChange"] > 1.58) #use abs() if consider all genes
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


###top domains for n4cel and cluster 3 genes in expressed cell types
stress3 = filter(stress,stress[,"log2FoldChange"] > 1.58)
index_dom = data.frame(c("GPCR|\\bSre\\b|Srg","\\bF-box|\\bFTH\\b|\\bHTH_48\\b|\\bFBA_2\\b","\\bHormone_recep|zf-C4\\b","\\bAcyl_transf_3|SGNH\\b","\\bBTB|BTB_2\\b"),c("GPCR","FBP","nhr","OAC","BTB"))
colnames(index_dom) = c("Var1","Freq")
type4 = c("n4cel","cls3")
exp_list2 = list(Cel_dup_N4Cel,cls3)
domain_stress = list()
for(i in 1:length(type4)){
  cur = stress3[stress3$WBGene %in% unique(exp_list2[[i]]$WBGene),]
  cur_dom = domain_elegans[domain_elegans$Wb %in% unique(cur$WBGene),]
  cur_summary = data.frame(table(cur_dom$domain))
  cur_summary2 = rbind(cur_summary,index_dom)
  for(j in 1:nrow(cur_summary2)){
    if(j <= (nrow(cur_summary2) -5)){
      cur_summary2[j,"domain_num"]= length(unique(cur_dom[grep(paste0("\\b",cur_summary2[j,"Var1"],"\\b"),cur_dom$domain),]$Wb))
      cur_summary2[j,"domain_total"]= length(unique(domain_elegans[grep(paste0("\\b",cur_summary2[j,"Var1"],"\\b"),domain_elegans$domain),]$Wb))
    } 
    else {
      cur_summary2[j,"domain_num"]= length(unique(cur_dom[grep(cur_summary2[j,"Var1"],cur_dom$domain),]$Wb))
      cur_summary2[j,"domain_total"]= length(unique(domain_elegans[grep(cur_summary2[j,"Var1"],domain_elegans$domain),]$Wb))
    }
  }
  cur_summary2$type = type4[i]
  domain_stress[[i]] = arrange(cur_summary2,desc(cur_summary2[,2]))  
}
domain_stress_2 = do.call(rbind,domain_stress)

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
  #cur = exp_list2[[i]]
  cur_L4 = L4_zscore[rownames(L4_zscore) %in% unique(cur$WBGene),]
  cur2 = filter(cur_L4,!is.na(cur_L4[,"AWB"]))
  cur2_mean = data.frame(sort(apply(cur2,2,mean2),decreasing = T))
  colnames(cur2_mean) = "L4"
  cur2_mean$type = type4[i]
  cur2_mean$celltype = rownames(cur2_mean)
  upregulate[[i]] = cur2_mean
}
upregulate2 = do.call(rbind,upregulate)

### L4 enrichment_N4cel
upregulate_N4cel = list()
for(i in 1:length(type4)){
  cur = exp_list2[[i]]
  cur_L4 = L4_zscore[rownames(L4_zscore) %in% unique(cur$WBGene),]
  cur2 = filter(cur_L4,!is.na(cur_L4[,"AWB"]))
  cur2_mean = data.frame(sort(apply(cur2,2,mean2),decreasing = T))
  colnames(cur2_mean) = "L4"
  cur2_mean$type = type4[i]
  cur2_mean$celltype = rownames(cur2_mean)
  upregulate_N4cel[[i]] = cur2_mean
}
upregulate_N4cel2 = do.call(rbind,upregulate_N4cel)


