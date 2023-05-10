library(stats)
library(purrr)
library(dplyr)
load("bulk_dev_RNAseq.RData")
load("RNAseq_conditions_FC.RData")
load("Cel_dup_N4Cel.RData")
load("Single_copy_genes.RData")
load("genes_dNdS.RData")
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

pdf("dNdS_n4cel_others_4types.pdf")
plot(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(n4cel_null$WBGene),]$w),verticals = TRUE,do.points =F,ylab="Cumulative Fraction",
     xlab="dN/dS",xlim=c(0,1),main="",col="#e3342f",lwd=3,cex.lab=1.5,cex.axis=2)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(n4cel_low$WBGene),]$w),verticals = TRUE, do.points =F,col="#f6993f",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(n4cel_medium$WBGene),]$w),verticals = TRUE, do.points =F,col="#ffed4a",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(n4cel_high$WBGene),]$w),verticals = TRUE, do.points =F,col="#38c172",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(genes_null2$WBGene),]$w),verticals = TRUE, do.points =F,col="#3490dc",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(genes_low2$WBGene),]$w),verticals = TRUE, do.points =F,col="#6574cd",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(genes_medium2$WBGene),]$w),verticals = TRUE, do.points =F,col="#9561e2",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(genes_high2$WBGene),]$w),verticals = TRUE, do.points =F,col="#f66d9b",lwd=3)
legend("bottomright",legend=c("other_genes_100","other_genes_10_100","other_genes_0_10","other_genes_0","n4cel_high","n4cel_medium","n4cel_low","n4cel_null"),lty=c(1,1),
       col=c("#f66d9b","#9561e2","#6574cd","#3490dc","#38c172","#ffed4a","#f6993f","#e3342f"),ncol=1,bty="n",inset=0.03,cex=1.6,lwd=4,text.font = 2)
dev.off()

pdf("dNdS_all_cls1_2_3_n4cel_single_copy.pdf")
plot(ecdf(genes_dNdS$w),verticals = TRUE,do.points =F,ylab="Cumulative Fraction",
     xlab="dN/dS",xlim=c(0,1),main="",col="#212121",lwd=3,cex.lab=1.5,cex.axis=2)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% rownames(cls1),]$w),verticals = TRUE, do.points =F,col="#D53E4F",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% rownames(cls3),]$w),verticals = TRUE, do.points =F,col="#3288BD",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% rownames(cls2),]$w),verticals = TRUE, do.points =F,col="#FFA138",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(Cel_dup_N4Cel$WBGene),]$w),verticals = TRUE, do.points =F,col="#5E4FA2",lwd=3)
lines(ecdf(genes_dNdS[genes_dNdS$wb %in% unique(genes_1to1$WBGene),]$w),verticals = TRUE, do.points =F,col="grey",lwd=3)
legend("bottomright",legend=c("All genes","cls1","cls3","cls2","n4cel_genes","single_copy_genes"),lty=c(1,1),col=c("#212121","#D53E4F","#3288BD","#FFA138","#5E4FA2","grey"),ncol=1,bty="n",inset=0.03,cex=1.6,lwd=6,text.font = 2)
dev.off()

#Fold change
RNAseq_conditions_FC3 = filter(RNAseq_conditions_FC,abs(RNAseq_conditions_FC[,"log2FoldChange"]) > 1.58)
NHR = domain_elegans[grep("\\bHormone_recep|zf-C4\\b",domain_elegans$domain),]
FBP = domain_elegans[grep("\\bF-box|\\bFTH\\b|\\bHTH_48\\b|\\bFBA_2\\b",domain_elegans$domain),]
MATH = domain_elegans[grep("\\bMATH\\b",domain_elegans$domain),]
LEC = domain_elegans[grep("\\bLectin_C\\b",domain_elegans$domain),]
domain_list = list(NHR,FBP,MATH,LEC)
domain_list2 = list()
type3 = c("NHR","FBP","MATH","LEC")
for(i in 1:length(type3)){
  cur = RNAseq_conditions_FC3[RNAseq_conditions_FC3$WBGene %in% unique(domain_list[[i]]$Wb),]
  cur$type = type3[i]
  domain_list2[[i]] = cur
}
domain_list2_2= do.call(rbind,domain_list2)
others = RNAseq_conditions_FC3[RNAseq_conditions_FC3$WBGene %in% setdiff(RNAseq_conditions_FC3$WBGene, domain_list2_2$WBGene),]
others$type = "others"
domain_list2_3 = rbind(domain_list2_2,others)
domain_list2_3$pvalue2 = -log10(domain_list2_3$pvalue)
write.csv(domain_list2_3,file = "domain_NHR_FBP_MATH_LEC_for_volcano.csv")


#enrichment of various types of genes with more than 3 foldchanges
RNAseq_conditions_FC2 = filter(RNAseq_conditions_FC,RNAseq_conditions_FC[,"log2FoldChange"] > 1.58) #use abs() if consider all genes
type = c("n4cel_null","n4cel_low","n4cel_medium","n4cel_high","cls1","cls2","cls3","other_genes_0","other_genes_0_10","other_genes_10_100","other_genes_100","single_copy_genes")
genes_1to1$type = "single-copy"
exp_list = list(n4cel_null,n4cel_low,n4cel_medium,n4cel_high,cls1,cls2,cls3,genes_null2,genes_low2,genes_medium2,genes_high2,genes_1to1)
exp_enr = data.frame(type)
for(i in 1:length(type)){
  cur = exp_list[[i]][exp_list[[i]]$WBGene %in% unique(RNAseq_conditions_FC2$WBGene),]
  exp_enr[i,"FC2"] = nrow(cur)
  exp_enr[i,"total"] = nrow(exp_list[[i]])
  
  exp_enr[i,"enrich"] = phyper(c(exp_enr[i,"FC2"]-1),exp_enr[i,"total"],c(19997 - exp_enr[i,"total"]),length(unique(RNAseq_conditions_FC2$WBGene)),lower.tail = F)
  exp_enr[i,"enrich2"] = -log10(exp_enr[i,"enrich"])
}
write.csv(exp_enr,file = "null_lowly_medium_high_FC3_enrichment_up_FC3.csv")

#plasticity
RNAseq_conditions_FC3 = select(RNAseq_conditions_FC,c("WBGene","log2FoldChange"))
RNAseq_conditions_FC3 = filter(RNAseq_conditions_FC3,log2FoldChange>0)
type = c("n4cel","n4cel_null","n4cel_low","n4cel_medium","n4cel_high","cls1","cls2","cls3","other_genes_0","other_genes_0_10","other_genes_10_100","other_genes_100","single_copy_genes")
exp_list = list(Cel_dup_N4Cel,n4cel_null,n4cel_low,n4cel_medium,n4cel_high,cls1,cls2,cls3,genes_null2,genes_low2,genes_medium2,genes_high2,genes_1to1)
exp_list2 = list()
for(i in 1:length(type)){
  cur = RNAseq_conditions_FC3[RNAseq_conditions_FC3$WBGene %in% unique(exp_list[[i]]$WBGene),]
  colnames(cur) = c("WBGene",type[i])
  exp_list2[[i]] = cur
}
exp_list2_3 = reduce(exp_list2,full_join,by = "WBGene")
write.csv(exp_list2_3,file = "FC_plasticity2_positive_FC.csv")

#versatility
RNAseq_conditions_FC3 = filter(RNAseq_conditions_FC,RNAseq_conditions_FC[,"log2FoldChange"] > 1.58)
type = c("n4cel","n4cel_null","n4cel_low","n4cel_medium","n4cel_high","cls1","cls2","cls3","other_genes_0","other_genes_0_10","other_genes_10_100","other_genes_100","single_copy_genes")
exp_list = list(Cel_dup_N4Cel,n4cel_null,n4cel_low,n4cel_medium,n4cel_high,cls1,cls2,cls3,genes_null2,genes_low2,genes_medium2,genes_high2,genes_1to1)
exp_list_versat = list()
exp_list_versat_mean = c()
for(i in 1:length(type)){
  cur = RNAseq_conditions_FC3[RNAseq_conditions_FC3$WBGene %in% unique(exp_list[[i]]$WBGene),]
  cur_sum = data.frame(table(cur$WBGene))
  colnames(cur_sum) = c("WBGene","Freq")
  cur_sum = filter(cur_sum,Freq >0)
  cur_sum$type = type[i]
  exp_list_versat[[i]] = cur_sum
  
  exp_list_versat_mean = append(exp_list_versat_mean,mean(cur_sum$Freq))
}
exp_list_versat_2 = do.call(rbind,exp_list_versat)
write.csv(exp_list_versat_2,file = "exp_list_versat.csv")

#the cell types of up-regulated genes
RNAseq_conditions_FC3 = filter(RNAseq_conditions_FC,abs(RNAseq_conditions_FC[,"log2FoldChange"]) > 1.58)
index_dom = data.frame(c("GPCR|\\bSre\\b|Srg","\\bF-box|\\bFTH\\b|\\bHTH_48\\b|\\bFBA_2\\b","\\bHormone_recep|zf-C4\\b","\\bAcyl_transf_3|SGNH\\b","\\bBTB|BTB_2\\b"),c("GPCR","FBP","nhr","OAC","BTB"))
colnames(index_dom) = c("Var1","Freq")
type4 = c("n4cel","cls2")
exp_list2 = list(Cel_dup_N4Cel,cls2)

domain_stress = list()
for(i in 1:length(type4)){
  cur = RNAseq_conditions_FC3[RNAseq_conditions_FC3$WBGene %in% unique(exp_list2[[i]]$WBGene),]
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
write.csv(domain_stress_2,file = "domains_upregulated_n4cel_cls2_FCmore3.csv")

### tissue enrichment
type2 = sort(colnames(L4_exp))
aaa2 = celltype[celltype$celltype2 %in% type2,]
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
  cur = RNAseq_conditions_FC3[RNAseq_conditions_FC3$WBGene %in% unique(exp_list2[[i]]$WBGene),]
  cur_L4 = L4_zscore[rownames(L4_zscore) %in% unique(cur$WBGene),]
  cur2 = filter(cur_L4,!is.na(cur_L4[,"AWB"]))
  cur2_mean = data.frame(sort(apply(cur2,2,mean2),decreasing = T))
  colnames(cur2_mean) = "L4"
  cur2_mean$type = type4[i]
  cur2_mean$celltype = rownames(cur2_mean)
  upregulate[[i]] = cur2_mean
}
upregulate2 = do.call(rbind,upregulate)
write.csv(upregulate2, file = "upregulate_celltypes_FC3.csv",row.names = F)

