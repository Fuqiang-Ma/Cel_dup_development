library(dplyr)
library(purrr)
load("cell_exp_domain.RData")
load("cluster_1_2_3.RData")
load("domain_elegans.RData")
load("celltype.RData")
load("neuronal_non neuronal_cells.RData")

#get neu_list and non_neu_list from "Neurons_non neuronal_cells_heatmap.R"
mean2 = function(x){
  mean(x,na.rm = T)
}

sd2 = function(x){
  sd(x,na.rm = T)
}

allcell_list = c(neu_list,non_neu_list)
type_allcell = c(type_neu,type_non_neu)
rec_z_cls2$WBGene = rownames(rec_z_cls2)
cls2 = filter(rec_z_cls2,cls ==2)
cls3 = filter(rec_z_cls2,cls ==3)

cls2_list = list()
for(i in 1:length(allcell_list)){
  cur = allcell_list[[i]][rownames(allcell_list[[i]]) %in% unique(cls2$WBGene),]
  cur$WBGene = rownames(cur)
  cls2_list[[i]] = cur
}

cls3_list = list()
for(i in 1:length(allcell_list)){
  cur = allcell_list[[i]][rownames(allcell_list[[i]]) %in% unique(cls3$WBGene),]
  cur$WBGene = rownames(cur)
  cls3_list[[i]] = cur
}

cls2_2 = filter(zscore_all,cls_type == "cls3" & X455 > 0.3)
cls3_2 = filter(zscore_all,cls_type == "cls2" & X3000 > 1)
cls2_3 = list(filter(cls2_2, celltype %in% c("BWM_head_row_B_VL")),filter(cls2_2, celltype %in% c("ASEL")))
cls3_3 = list(filter(cls3_2, celltype %in% c("hyp1_2")),filter(cls3_2, celltype %in% c("ASG")))

#cluster 2 genes
cells171$time_range = paste("X_",cells171$Start_time,cells171$End_time,sep = "_")
cls2_list2 = do.call(rbind,cls2_list)
type = c("BWM_head_row_B_VL","ASEL")
for(i in c(1:2)){
  cur = cls2_list2[cls2_list2$celltype %in% unique(cls2_3[[i]]$celltype),]
  cur_domain = domain_elegans[domain_elegans$Wb %in% unique(cur$WBGene),]
  cur_domain_sum = data.frame(table(cur_domain$domain))
  index_dom = data.frame(c("GPCR|\\bSre\\b|Srg","\\bF-box|\\bFTH\\b|\\bHTH_48\\b|\\bFBA_2\\b","\\bHormone_recep\\b|\\bzf-C4\\b","\\bAcyl_transf_3\\b|\\bSGNH\\b","\\bBTB\\b|\\bBTB_2\\b"),c("GPCR","FBP","nhr","OAC","BTB"))
  colnames(index_dom) = c("Var1","Freq")
  cur_domain_sum2 = rbind(cur_domain_sum,index_dom)
  for( m in 1:nrow(cur_domain_sum2)){
    if(m <= (nrow(cur_domain_sum2) -5)){
      cur_domain_sum2[m,"domain_num"]= length(unique(cur_domain[grep(paste0("\\b",cur_domain_sum2[m,"Var1"],"\\b"),cur_domain$domain),]$Wb))
      cur_domain_sum2[m,"domain_total"]= length(unique(domain_elegans[grep(paste0("\\b",cur_domain_sum2[m,"Var1"],"\\b"),domain_elegans$domain),]$Wb))
    } 
    else {
      cur_domain_sum2[m,"domain_num"]= length(unique(cur_domain[grep(cur_domain_sum2[m,"Var1"],cur_domain$domain),]$Wb))
      cur_domain_sum2[m,"domain_total"]= length(unique(domain_elegans[grep(cur_domain_sum2[m,"Var1"],domain_elegans$domain),]$Wb))
    }
  }
  cur_domain_sum3 = arrange(cur_domain_sum2,desc(domain_num))
  cls2_list3 = list()
  for(n in 1:10){
    dom = domain_elegans[grep(cur_domain_sum3[n,"Var1"],domain_elegans$domain),]
    dom2 = cur[cur$WBGene %in% unique(dom$Wb),]
    dom2_mean = data.frame(apply(dom2[,1:77],2,mean2),apply(dom2[,1:77],2,sd2))
    colnames(dom2_mean) = c(paste0(cur_domain_sum3[n,"Var1"],"_mean(",cur_domain_sum3[n,"domain_num"],")"),paste0(cur_domain_sum3[n,"Var1"],"_sd"))
    dom2_mean = arrange(dom2_mean,as.numeric(rownames(dom2_mean)))
    cls2_list3[[n]] = dom2_mean
  }
  cls2_list3_2 = do.call(bind_cols,cls2_list3)
  write.csv(cls2_list3_2,file = paste0("cls2_",type[i],"_domain.csv"))
  
  #lowly expressed
  type2 = c("hyp1_2","ASG")
  cls2_2_less0 = filter(zscore_all,cls_type == "cls3" & X455 < 0)
  cls2_3_less0 = list(filter(cls2_2_less0, celltype %in% c("hyp1_2")),filter(cls2_2_less0, celltype %in% c("ASG")))
  cur2 = allcell_list[which(type_allcell %in% as.character(unique(cls2_3_less0[[i]]$celltype)))]
  for(o in 1:length(cur2)){
    cur2[[o]]$WBGene = rownames(cur2[[o]])
  }
  cur_2 = do.call(rbind,cur2)
  cls2_list3_3 = list()
  for(p in 1:10){
    dom_low = domain_elegans[grep(cur_domain_sum3[p,"Var1"],domain_elegans$domain),]
    dom3 = cur[cur$WBGene %in% unique(dom_low$Wb),]
    dom2_low = cur_2[cur_2$WBGene %in% unique(dom3$WBGene),]
    dom2_mean_low = data.frame(apply(dom2_low[,1:77],2,mean2),apply(dom2_low[,1:77],2,sd2))
    colnames(dom2_mean_low) = c(paste0(cur_domain_sum3[p,"Var1"],"_mean(",cur_domain_sum3[p,"domain_num"],")"),paste0(cur_domain_sum3[p,"Var1"],"_sd"))
    dom2_mean_low = arrange(dom2_mean_low,as.numeric(rownames(dom2_mean_low)))
    cls2_list3_3[[p]] = dom2_mean_low
  }
  cls2_list3_4 = do.call(bind_cols,cls2_list3_3)
  write.csv(cls2_list3_4,file = paste0("cls2_lowexp_",type2[i],"_domain.csv"))
}


#cluster 3 genes
cls3_list2 = do.call(rbind,cls3_list)
for(i in c(1:2)){
  type = c("hyp1_2","ASG")
  cur = cls3_list2[cls3_list2$celltype %in% unique(cls3_3[[i]]$celltype),]
  cur_domain = domain_elegans[domain_elegans$Wb %in% unique(cur$WBGene),]
  cur_domain_sum = data.frame(table(cur_domain$domain))
  index_dom = data.frame(c("GPCR|\\bSre\\b|Srg","\\bF-box|\\bFTH\\b|\\bHTH_48\\b|\\bFBA_2\\b","\\bHormone_recep\\b|\\bzf-C4\\b","\\bAcyl_transf_3\\b|\\bSGNH\\b","\\bBTB\\b|\\bBTB_2\\b"),c("GPCR","FBP","nhr","OAC","BTB"))
  colnames(index_dom) = c("Var1","Freq")
  cur_domain_sum2 = rbind(cur_domain_sum,index_dom)
  for( m in 1:nrow(cur_domain_sum2)){
    if(m <= (nrow(cur_domain_sum2) -5)){
      cur_domain_sum2[m,"domain_num"]= length(unique(cur_domain[grep(paste0("\\b",cur_domain_sum2[m,"Var1"],"\\b"),cur_domain$domain),]$Wb))
      cur_domain_sum2[m,"domain_total"]= length(unique(domain_elegans[grep(paste0("\\b",cur_domain_sum2[m,"Var1"],"\\b"),domain_elegans$domain),]$Wb))
    } 
    else {
      cur_domain_sum2[m,"domain_num"]= length(unique(cur_domain[grep(cur_domain_sum2[m,"Var1"],cur_domain$domain),]$Wb))
      cur_domain_sum2[m,"domain_total"]= length(unique(domain_elegans[grep(cur_domain_sum2[m,"Var1"],domain_elegans$domain),]$Wb))
    }
  }
  cur_domain_sum3 = arrange(cur_domain_sum2,desc(domain_num))
  cls3_list3 = list()
  for(n in 1:10){
    dom = domain_elegans[grep(cur_domain_sum3[n,"Var1"],domain_elegans$domain),]
    dom2 = cur[cur$WBGene %in% unique(dom$Wb),]
    dom2_mean = data.frame(apply(dom2[,1:77],2,mean2),apply(dom2[,1:77],2,sd2))
    colnames(dom2_mean) = c(paste0(cur_domain_sum3[n,"Var1"],"_mean(",cur_domain_sum3[n,"domain_num"],")"),paste0(cur_domain_sum3[n,"Var1"],"_sd"))
    dom2_mean = arrange(dom2_mean,as.numeric(rownames(dom2_mean)))
    cls3_list3[[n]] = dom2_mean
  }
  cls3_list3_2 = do.call(bind_cols,cls3_list3)
  write.csv(cls3_list3_2,file = paste0("cls3_",type[i],"_domain.csv"))
  
  #lowly expressed
  type2 = c("BWM_head_row_B_VL","ASEL")
  cls3_2_less0 = filter(zscore_all,cls_type == "cls2" & X3000 < 0.6)
  cls3_3_less0 = list(filter(cls3_2_less0, celltype %in% c("BWM_head_row_B_VL")),filter(cls3_2_less0, celltype %in% c("ASEL")))
  cur2 = allcell_list[which(type_allcell %in% as.character(unique(cls3_3_less0[[i]]$celltype)))]
    for(o in 1:length(cur2)){
      cur2[[o]]$WBGene = rownames(cur2[[o]])
    }
    cur_2 = do.call(rbind,cur2)
    cls3_list3_3 = list()
    for(p in 1:10){
      dom_low = domain_elegans[grep(cur_domain_sum3[p,"Var1"],domain_elegans$domain),]
      dom3 = cur[cur$WBGene %in% unique(dom_low$Wb),]
      dom2_low = cur_2[cur_2$WBGene %in% unique(dom3$WBGene),]
      dom2_mean_low = data.frame(apply(dom2_low[,1:77],2,mean2),apply(dom2_low[,1:77],2,sd2))
      colnames(dom2_mean_low) = c(paste0(cur_domain_sum3[p,"Var1"],"_mean(",cur_domain_sum3[p,"domain_num"],")"),paste0(cur_domain_sum3[p,"Var1"],"_sd"))
      dom2_mean_low = arrange(dom2_mean_low,as.numeric(rownames(dom2_mean_low)))
      cls3_list3_3[[p]] = dom2_mean_low
    }
    cls3_list3_4 = do.call(bind_cols,cls3_list3_3)
    write.csv(cls3_list3_4,file = paste0("cls3_lowexp_",type2[i],"_domain.csv"))
}

