library(dplyr)
library(stringr)
load("sc_exp.RData")
load("celltype.RData")
load("cluster_1_2_3.RData")
load("cell_type_num.RData")

### all lineage cells
tissue_lin = list()
tissue = c("AB","MS","D","C","E","P4")
for(i in 1:length(tissue)){
  lin_2 = filter(cell_num_lineage, lin2 == tissue[i])
  lin_type = sort(unique(lin_2$Cell))
  lin_list=list()
  for(j in 1:length(lin_type)){
    cur = filter(lin_2,Cell == lin_type[j])
    cur_type = celltype[celltype$celltype2 %in% cur$Annotated.lineage_Packer,]
    cur_exp = select(sc_exp,c(cur_type[,"index"])) 
    
    for(g in 1:ncol(cur_exp)){
      cur2 = filter(cur_type, index == colnames(cur_exp)[g])
      colnames(cur_exp)[g] = cur2$celltype2
      
      cur3 = filter(cur,Annotated.lineage_Packer == colnames(cur_exp)[g])
      colnames(cur_exp)[g] = cur3$time_range
    }
    lin_list[[j]] = cur_exp
  }
  
  time_lin = c(seq(55,510,by = 10))
  lin_list2 = list()
  for(o in 1:length(lin_list)){
    df = data.frame(matrix(nrow = nrow(lin_list[[o]]), ncol = length(time_lin) + 1))
    rownames(df) = rownames(lin_list[[o]])
    colnames(df) = c(time_lin,"celltype")
    df[,"celltype"] = lin_type[o] 
    for(p in 1:ncol(lin_list[[o]])){
      for(m in 1:(ncol(df) -1)){
        if(isTRUE(as.numeric(colnames(df)[m]) > as.numeric(strsplit(colnames(lin_list[[o]])[p],"_")[[1]][2]) & as.numeric(colnames(df)[m]) < as.numeric(strsplit(colnames(lin_list[[o]])[p],"_")[[1]][3]))){
          df[,m] = lin_list[[o]][,p]
        }
      }
    }
    lin_list2[[o]] = df
  }
  
  for(q in 1:length(lin_list2)){
    lin_list2[[q]]$WBGene = rownames(lin_list2[[q]])
  }
  
  lin_list2_2 = do.call(rbind,lin_list2)
  lin_all = data.frame(matrix(nrow = nrow(lin_list2[[1]]), ncol = length(time_lin)))
  rownames(lin_all) = rownames(lin_list2[[1]])
  colnames(lin_all) = c(time_lin)
  aaaa = lin_2[!duplicated(lin_2$Cell),]
  sum(as.numeric(aaaa$X))
  for(s in 1:nrow(lin_all)){
    cur = filter(lin_list2_2,WBGene == rownames(lin_all)[s])
    cur2 = select(cur,-c("celltype","WBGene"))
    lin_all[s,] = apply(cur2,2,function(x){ifelse(all(is.na(x)),NA,sum(!is.na(x))*mean(x,na.rm=TRUE))})
    if(s %% 10000 ==0){
      cat(s,"rows finished \n")
    }
  }
  lin_num = c()
  for(t in 1:length(time_lin)){
    cur= filter(lin_2,time1 < time_lin[t] & time2 > time_lin[t])
    lin_num = append(lin_num,sum(as.numeric(cur$X)))
  }
  for(u in 1:length(lin_num)){
    if(isTRUE(lin_num[u] == 0)){
      lin_num[u] =1
    }
  }
  for(v in 1:ncol(lin_all)){
    lin_all[,v] = lin_all[,v]/lin_num[v]
  }
  tissue_lin[[i]] = lin_all
}

#late embryonic cells
tissue_emb = list()
tissue = c("AB","MS","D","C","E","P4")
for(i in 1:length(tissue)){
  emb_late3_2 = filter(cell_num_emb, lin2 == tissue[i])
  emb_type = sort(unique(emb_late3_2$type))
  emb_list=list()
  for(j in 1:length(emb_type)){
    cur = filter(emb_late3_2,type == emb_type[j])
    cur_type = celltype[celltype$celltype2 %in% cur$type2,]
    cur_exp = select(sc_exp,c(cur_type[,"index"])) 
    for(g in 1:ncol(cur_exp)){
      cur2 = filter(cur_type, index == colnames(cur_exp)[g])
      colnames(cur_exp)[g] = cur2$celltype2
      
      cur3 = filter(emb_late3_2,type2 == colnames(cur_exp)[g])
      cur_exp[,g] = data.frame(as.numeric(cur_exp[,g])) * as.numeric(cur3$X)
    }
    emb_list[[j]] = cur_exp
  }

  emb_late3_2$time_range = paste("X",emb_late3_2$time1,emb_late3_2$time2,sep = "_")
  for(h in 1:length(emb_list)){
    for(k in 1:ncol(emb_list[[h]])){
      for(n in 1:nrow(emb_late3_2)){
        if(colnames(emb_list[[h]])[k] == emb_late3_2[n,"type2"]){
          colnames(emb_list[[h]])[k] = emb_late3_2[n,"time_range"]
        }
      }
    }
  }
  
  time_emb = c(seq(515,800,by = 10))
  emb_list2 = list()
  for(o in 1:length(emb_list)){
    df = data.frame(matrix(nrow = nrow(emb_list[[o]]), ncol = length(time_emb) + 1))
    rownames(df) = rownames(emb_list[[o]])
    colnames(df) = c(time_emb,"celltype")
    df[,"celltype"] = emb_type[o] 
    for(p in 1:ncol(emb_list[[o]])){
      for(m in 1:(ncol(df) -1)){
        if(isTRUE(as.numeric(colnames(df)[m]) > as.numeric(strsplit(colnames(emb_list[[o]])[p],"_")[[1]][2]) & as.numeric(colnames(df)[m]) < as.numeric(strsplit(colnames(emb_list[[o]])[p],"_")[[1]][3]))){
          df[,m] = emb_list[[o]][,p]
        }
      }
    }
    emb_list2[[o]] = df
  }
  
  for(q in 1:length(emb_list2)){
    emb_list2[[q]]$WBGene = rownames(emb_list2[[q]])
  }
  
  emb_list2_2 = do.call(rbind,emb_list2)
  emb_all = data.frame(matrix(nrow = nrow(emb_list2[[1]]), ncol = length(time_emb)))
  rownames(emb_all) = rownames(emb_list2[[1]])
  colnames(emb_all) = c(time_emb)
  aaaa = emb_late3_2[!duplicated(emb_late3_2$type),]
  sum(as.numeric(aaaa$X))
  for(s in 1:nrow(emb_all)){
    cur = filter(emb_list2_2,WBGene == rownames(emb_all)[s])
    cur2 = select(cur,-c("celltype","WBGene"))
    emb_all[s,] = apply(cur2,2,function(x){ifelse(all(is.na(x)),NA,sum(!is.na(x))*mean(x,na.rm=TRUE))}) #can also use apply(cur2,2,sum,na.rm=T),but will give NA if all valuses are NA!
    if(s %% 10000 ==0){
      cat(s,"rows finished \n")
    }
  }
  emb_num = c()
  for(t in 1:length(time_emb)){
    cur= filter(emb_late3_2,time1 < time_emb[t] & time2 > time_emb[t])
    emb_num = append(emb_num,sum(as.numeric(cur$X)))
  }
  for(u in 1:length(emb_num)){
    if(isTRUE(emb_num[u] == 0)){
      emb_num[u] =1
    }
  }
  for(v in 1:ncol(emb_all)){
    emb_all[,v] = emb_all[,v]/emb_num[v] # it will produce error if use below code!!
  }
  tissue_emb[[i]] = emb_all 
}


###L2 cells
tissue_l2 = list()
tissue = c("AB","MS","D","C","E","P4")
for(i in 1:length(tissue)){
  l2_2 = filter(cell_num_L2, lin2 == tissue[i])
  l2_type = sort(unique(l2_2$X))
  l2_list=list()
  for(j in 1:length(l2_type)){
    cur = filter(l2_2,X == l2_type[j])
    cur_type = celltype[celltype$celltype2 %in% cur$type,]
    cur_exp = select(sc_exp,c(cur_type[,"index"]))
    for(g in 1:ncol(cur_exp)){
      cur2 = filter(cur_type, index == colnames(cur_exp)[g])
      colnames(cur_exp)[g] = cur2$celltype2
      
      cur3 = filter(l2_2,type == colnames(cur_exp)[g])
      cur_exp[,g] = data.frame(as.numeric(cur_exp[,g])) * as.numeric(cur3$X.1)
    }
    l2_list[[j]] = cur_exp
  }
  
  for(q in 1:length(l2_list)){
    colnames(l2_list[[q]]) = "L2"
    l2_list[[q]]$WBGene = rownames(l2_list[[q]])
  }
  
  l2_list_2 = do.call(rbind,l2_list)
  l2_all = data.frame(matrix(nrow = nrow(l2_list[[1]]), ncol = 1))
  rownames(l2_all) = rownames(l2_list[[1]])
  colnames(l2_all) = "L2"
  aaaa = l2_2[!duplicated(l2_2$type),]
  sum(as.numeric(aaaa$X.1))
  for(s in 1:nrow(l2_all)){
    cur = filter(l2_list_2,WBGene == rownames(l2_all)[s])
    cur2 = select(cur,-c("WBGene"))
    l2_all[s,] = apply(cur2,2,function(x){ifelse(all(is.na(x)),NA,sum(!is.na(x))*mean(x,na.rm=TRUE))})
    if(s %% 10000 ==0){
      cat(s,"rows finished \n")
    }
  }
  l2_all2 = l2_all/sum(as.numeric(aaaa$X.1))
  tissue_l2[[i]] = l2_all2
}


###L4 cells
tissue_l4 = list()
tissue = c("AB","MS","D","C","E","P4")
for(i in 1:length(tissue)){
  l4_2 = filter(cell_num_L4, lin2 == tissue[i])
  l4_type = sort(unique(l4_2$celltype))
  l4_list=list()
  for(j in 1:length(l4_type)){
    cur = filter(l4_2,celltype == l4_type[j])
    cur_type = celltype[celltype$celltype2 %in% cur$celltype,]
    cur_exp = select(sc_exp,c(cur_type[,"index"]))
    for(g in 1:ncol(cur_exp)){
      cur2 = filter(cur_type, index == colnames(cur_exp)[g])
      colnames(cur_exp)[g] = cur2$celltype2
      
      cur3 = filter(l4_2,celltype == colnames(cur_exp)[g])
      cur_exp[,g] = data.frame(as.numeric(cur_exp[,g])) * as.numeric(cur3$X)
    }
    l4_list[[j]] = cur_exp
  }
  
  for(q in 1:length(l4_list)){
    colnames(l4_list[[q]]) = "L4"
    l4_list[[q]]$WBGene = rownames(l4_list[[q]])
  }
  
  l4_list_2 = do.call(rbind,l4_list)
  l4_all = data.frame(matrix(nrow = nrow(l4_list[[1]]), ncol = 1))
  rownames(l4_all) = rownames(l4_list[[1]])
  colnames(l4_all) = "L4"
  aaaa = l4_2[!duplicated(l4_2$celltype),]
  sum(as.numeric(aaaa$X))
  for(s in 1:nrow(l4_all)){
    cur = filter(l4_list_2,WBGene == rownames(l4_all)[s])
    cur2 = select(cur,-c("WBGene"))
    l4_all[s,] = apply(cur2,2,function(x){ifelse(all(is.na(x)),NA,sum(!is.na(x))*mean(x,na.rm=TRUE))})
    if(s %% 10000 ==0){
      cat(s,"rows finished \n")
    }
  }
  l4_all2 = l4_all/sum(as.numeric(aaaa$X))
  tissue_l4[[i]] = l4_all2
}

#exp
all_list = c(tissue_lin,tissue_emb,tissue_l2,tissue_l4)
cell_list= list()
for(i in 1:6){
  num = seq(i,25,by=6)
  cur = all_list[num]
  cur = cur[!sapply(cur, is.null)]
  for(j in 1:length(cur)){
    cur[[j]] = arrange(cur[[j]],rownames(cur[[j]]))
  }
  cur2 = do.call(bind_cols,cur)
  cell_list[[i]] = cur2
}

tissue = c("AB","MS","D","C","E","P4")
cls_list = list()
for(i in 1:3){
  cur = filter(rec_z_cls2, cls == i)
  cls_list2 = list()
  for(j in 1:length(tissue)){
    cur2 = log2(cell_list[[j]]+1)
    cur_cell = cur2[rownames(cur),]
    cur3 = apply(cur_cell,2,mean2)
    cur3 = data.frame(t(cur3))
    cur3$type = tissue[j]
    cls_list2[[j]] = cur3
  }
  cls_list2_2 = data.frame(do.call(rbind,cls_list2))
  cls_list2_2$type2 = paste0("cls_",i)
  cls_list[[i]] = cls_list2_2
}
cls_list_2 = do.call(rbind,cls_list)
cls_list_2 <- apply(cls_list_2,2,as.character)
write.csv(cls_list_2, file = "cells6_exp_cls_1_2_3.csv")

