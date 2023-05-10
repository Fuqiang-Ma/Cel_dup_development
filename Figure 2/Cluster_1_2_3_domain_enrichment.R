library(stats)
library(dplyr)
load("domain_elegans.RData")
load("cluster_1_2_3.RData")

index_dom = data.frame(c("GPCR|\\bSre\\b|Srg","\\bF-box|FTH|HTH_48|FBA_2\\b","\\bHormone_recep|zf-C4\\b","\\bAcyl_transf_3|SGNH\\b","\\bBTB|BTB_2\\b"),c("GPCR","FBP","nhr","OAC","BTB"))
colnames(index_dom) = c("Var1","Freq")
cls_domain = list()
for( i in 1:3){
  cls = filter(rec_z_cls2,cls == i) 
  cur = domain_elegans[domain_elegans$Wb %in% rownames(cls),]
  cur_summary = data.frame(table(cur$domain))
  cur_summary2 = rbind(cur_summary,index_dom)
  for(j in 1:nrow(cur_summary2)){
    if(j <= (nrow(cur_summary2) -5)){
      cur_summary2[j,"domain_num"]= length(unique(cur[grep(paste0("\\b",cur_summary2[j,"Var1"],"\\b"),cur$domain),]$Wb))
      cur_summary2[j,"domain_total"]= length(unique(domain_elegans[grep(paste0("\\b",cur_summary2[j,"Var1"],"\\b"),domain_elegans$domain),]$Wb))
    } 
    else {
      cur_summary2[j,"domain_num"]= length(unique(cur[grep(cur_summary2[j,"Var1"],cur$domain),]$Wb))
      cur_summary2[j,"domain_total"]= length(unique(domain_elegans[grep(cur_summary2[j,"Var1"],domain_elegans$domain),]$Wb))
    }
  }
  enr = function(x,y){
    z = phyper(x-1, nrow(cls), c(19997 - nrow(cls)), y, lower.tail=F) #example: phyper(169-1, 1000, 1183, 261, lower.tail=F)
  }
  for(m in 1:nrow(cur_summary2)){
    cur_summary2[m,"enrichment"] = enr(cur_summary2[m,"domain_num"],cur_summary2[m,"domain_total"])
    cur_summary2[m,"enrichment2"] = -log10(cur_summary2[m,"enrichment"])
  }
  colnames(cur_summary2) = paste0(colnames(cur_summary2),"_cls",i)
  cls_domain[[i]] = arrange(cur_summary2,desc(cur_summary2[,5]))
}
enrichment = merge(cls_domain[[1]],cls_domain[[2]],by.x = "Var1_cls1", by.y = "Var1_cls2",all = T)
enrichment2 = merge(enrichment,cls_domain[[3]],by.x = "Var1_cls1", by.y = "Var1_cls3", all = T)
write.csv(enrichment2, file = "enrichment_all_cluster1_2_3.csv",row.names = F)
