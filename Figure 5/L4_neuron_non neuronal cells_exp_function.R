library(tidyr)
library(dplyr)
load("L4_exp.RData")
load("domain_elegans.RData")
load("cluster_1_2_3.RData")
load("Cel_dup_N4Cel.RData")
load("Single_copy_genes.RData")
load("gene_count_orthologs.RData")

L4_exp2 = data.frame(t(L4_exp))
for(i in 1:ncol(L4_exp2)){
  L4_exp2[165,i] = nrow(filter(L4_exp2[1:164,],as.numeric(L4_exp2[1:164,i]) >0))
  if(i %% 3000 == 0){
    cat(i," rows finished \n")
  }
}

#gene exp in cell type
L4_exp2 = data.frame(t(L4_exp2))
from3 = seq(1,166,by=10)
to3 = append(seq(10,166,by=10),164)
type3 = c("n1_10","n11_20","n21_30","n31_40","n41_50","n51_60","n61_70","n71_80","n81_90","n91_100","n101_110","n111_120","n121_130","n131_140","n141_150","n151_160","n161_164")
index = data.frame(from3, to3, type3)
for ( i in 1:nrow(L4_exp2)){
  for(j in 1:nrow(index)){
    if(as.numeric(L4_exp2[i,"X165"]) >= index[j,"from3"] & as.numeric(L4_exp2[i,"X165"]) <= index[j,"to3"]){
      L4_exp2[i,"type"] = index[j,"type3"]
    }
  }
}

#enrichment
cls3 = filter(rec_z_cls2,cls == 3)
types_enr = data.frame(rev(sort(table(L4_exp2$type))))
types_enr2 = merge(types_enr,index,by.x = "Var1", by.y ="type3")
types_enr2 = arrange(types_enr2,from3)

for(i in 1:nrow(types_enr2)){
  cur = filter(L4_exp2,type == types_enr2[i,"Var1"])
  
  ###n4cel genes
  cur_dup = cur[unique(intersect(rownames(cur),Cel_dup_N4Cel$WBGene)),]
  types_enr2[i,"n4cel"] = nrow(cur_dup)
  enr = function(x){
    z = phyper(x-1, types_enr2[i,"Freq"], 19069 - types_enr2[i,"Freq"], 4483, lower.tail=F)
  }
  types_enr2[i,"enrichment"] = enr(nrow(cur_dup))
  types_enr2[i,"enrichment2"] = -log10(types_enr2[i,"enrichment"])
  
  ###single copy genes
  cur_single = cur[unique(intersect(rownames(cur),genes_1to1$WBGene)),]
  types_enr2[i,"single_copy"] = nrow(cur_single)
  enr_signle = function(x){
    z = phyper(x-1, types_enr2[i,"Freq"], 19069 - types_enr2[i,"Freq"], 4720, lower.tail=F)
  }
  types_enr2[i,"enrichment_single"] = enr_signle(nrow(cur_single))
  types_enr2[i,"enrichment2_single"] = -log10(types_enr2[i,"enrichment_single"])
  
  ###cls3 genes
  cur_cls3 = cur[unique(intersect(rownames(cur),rownames(cls3))),]
  types_enr2[i,"cls3"] = nrow(cur_cls3)
  enr_signle = function(x){
    z = phyper(x-1, types_enr2[i,"Freq"], 19069 - types_enr2[i,"Freq"], 781, lower.tail=F)
  }
  types_enr2[i,"enrichment_cls3"] = enr_signle(nrow(cur_cls3))
  types_enr2[i,"enrichment_cls3"] = -log10(types_enr2[i,"enrichment_cls3"])
}
write.csv(types_enr2, file = "enrichment_single_copy_cls3_n4cel_genes.csv")


#enrichment for recently duplicated genes
index_dom = data.frame(c("GPCR|\\bSre\\b|Srg","\\bF-box|FTH|HTH_48|FBA_2\\b","\\bHormone_recep|zf-C4\\b","\\bAcyl_transf_3|SGNH\\b","\\bBTB|BTB_2\\b"),c("GPCR","FBP","nhr","OAC","BTB"))
colnames(index_dom) = c("Var1","Freq")
types_domain_enr = data.frame(rev(sort(table(L4_exp2$type))))
types_domain_enr2 = merge(types_domain_enr,index,by.x = "Var1", by.y ="type3")
types_domain_enr2 = arrange(types_domain_enr2,from3)

library(stats)
domain_enr = list()
for(i in 1:nrow(types_domain_enr2)){
  cur = filter(L4_exp2, type == types_domain_enr2[i,"Var1"])
  cur = cur[unique(intersect(rownames(cur),Cel_dup_N4Cel$WBGene)),]
  cur_dom = domain_elegans[domain_elegans$Wb %in% rownames(cur),]
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
  enr = function(x,y){
    z = phyper(x-1, types_domain_enr2[i,"Freq"], c(19997 - types_domain_enr2[i,"Freq"]), y, lower.tail=F)
  }
  for(m in 1:nrow(cur_summary2)){
    cur_summary2[m,"enrichment"] = enr(cur_summary2[m,"domain_num"],cur_summary2[m,"domain_total"])
    cur_summary2[m,"enrichment2"] = -log10(cur_summary2[m,"enrichment"])
  }
  cur_summary2$type = types_domain_enr2[i,"Var1"]
  domain_enr[[i]] = arrange(cur_summary2,desc(cur_summary2[,6]))
}
domain_enr_n4cel = do.call(rbind,domain_enr)
write.csv(domain_enr_n4cel,file = "domain_enr_n4cel.csv")

#enrichment for signle copy genes
domain_enr_single = list()
aaa3 = data.frame("n141_164",sum(as.numeric(types_domain_enr2[15:17,"Freq"])),"111","164")
colnames(aaa3) = colnames(types_domain_enr2)
for(i in 1:nrow(aaa3)){
  cur = L4_exp2[L4_exp2$type %in% c("n141_150","n151_160","n161_164"),]
  cur = cur[unique(intersect(rownames(cur),genes_1to1$WBGene)),]
  cur_dom = domain_elegans[domain_elegans$Wb %in% rownames(cur),]
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
  enr = function(x,y){
    z = phyper(x-1, aaa3[i,"Freq"], c(19997 - aaa3[i,"Freq"]), y, lower.tail=F)
  }
  for(m in 1:nrow(cur_summary2)){
    cur_summary2[m,"enrichment"] = enr(cur_summary2[m,"domain_num"],cur_summary2[m,"domain_total"])
    cur_summary2[m,"enrichment2"] = -log10(cur_summary2[m,"enrichment"])
  }
  cur_summary2$type = aaa3[i,"Var1"]
  domain_enr_single[[i]] = arrange(cur_summary2,desc(cur_summary2[,6]))  
}
domain_enr_single_copy = do.call(rbind,domain_enr_single)
write.csv(domain_enr_single_copy,file = "domain_enr_single_copy_n141_164.csv")

# number of genes in specific neurons and non-neuronal cells
n1_10 = filter(L4_exp2,type == "n1_10")
for(i in 1:(ncol(L4_exp2) -2)){
  cur = filter(L4_exp2,as.numeric(L4_exp2[,i]) > 0)
  cur_n4cel = cur[unique(intersect(rownames(cur),duplication_N4Cel399$Gene)),]
  
  n1_10_n4cel = length(unique(intersect(rownames(cur_n4cel),rownames(n1_10))))
  n1_10_all = length(unique(intersect(rownames(n1_10),rownames(cur))))
  
  L4_exp2[19070,i] = n1_10_n4cel
  L4_exp2[19071,i] = nrow(cur_n4cel)
  L4_exp2[19072,i] = nrow(cur)
  L4_exp2[19073,i] = n1_10_all
}
n4cel_n1_10 = data.frame(t(L4_exp2[19070:19073,1:164]))
colnames(n4cel_n1_10) = c("#_n1_10_n4cel", "#_exp_n4cel","#_total_exp_genes","#_n1_10_genes")
write.csv(n4cel_n1_10,file = "num_genes_each_cell_n4cel_n1_10.csv")


# num of ortholog genes
type4 = c("n1_10","n11_20","n21_30","n31_40","n41_50","n51_60","n61_70","n71_80","n81_90","n91_100","n101_110","n111_120","n121_130","n131_140","n141_150","n151_160","n161_164")
plot = list()
for(i in 1:length(type4)){   
  cur = filter(L4_exp2,type == type4[i])
  cur2 = gene_count_orthologs[gene_count_orthologs$WBGene %in% rownames(cur),]
  plot[[i]] = cur2
}

species = c("C. becei","C. bovis","C. briggsae","C. elegans","C. inopinata","C. latens","C. nigoni","C. panamensis","C. remanei","C. tribulationis","C. tropicalis")
species2 =  list()
for(i in 1:length(species)){
  ortho_num = data.frame(matrix(nrow = 4,ncol = length(type4)))
  rownames(ortho_num)=c("0","1","2",">2")
  colnames(ortho_num)= type4
  for(j in 1:length(plot)){    
    ortho_num[1,j] = nrow(filter(plot[[j]],plot[[j]][,species[i]] ==0))
    ortho_num[2,j] = nrow(filter(plot[[j]],plot[[j]][,species[i]] ==1))
    ortho_num[3,j] = nrow(filter(plot[[j]],plot[[j]][,species[i]] ==2))
    ortho_num[4,j] = nrow(filter(plot[[j]],plot[[j]][,species[i]] > 2))
  }
  ortho_num2 = data.frame(t(ortho_num))
  colnames(ortho_num2)=c("0","1","2",">2")
  ortho_num2[,"types14"] = rownames(ortho_num2)
  ortho_num2[,"species_type"] = species[i]
  species2[[i]] =  ortho_num2
}
species2_2 = do.call(rbind,species2)
write.csv(species2_2,file="spceis11_types14_homolog_num.csv",row.names=F)

