library(dplyr)
library(RColorBrewer)
library(gplots)
library(stringr)
library(purrr)
load("cluster_1_2_3.RData")
load("sc_exp.RData")
load("celltype.RData")
load("neuronal_non neuronal_cells.RData")

zscore = function(x){
  x = (x - mean(x))/sd(x)
}

mean2 = function(x){
  mean(x,na.rm = T)
}


#neurons
neu_list2 = list()
for(i in 1:length(type_neu)){
  for(j in 1:nrow(celltype)){
    cur = filter(neu_cell,Neuron == type_neu[i])
    cur_type = celltype[celltype$celltype2 %in% cur$Precursor2,]
    cur_exp = select(sc_exp,c(cur_type[,"index"])) #note whether the order is time series
    for(g in 1:ncol(cur_exp)){
      for(h in 1:nrow(cur_type)){
        if(colnames(cur_exp)[g] == cur_type[h,"index"]){
          colnames(cur_exp)[g] = cur_type[h,"celltype2"]
        }}}
  }
  neu_list2[[i]] = cur_exp
  if(i %% 10 ==0){
    cat(i,"rows finished \n")
  }
}
neu_list2[[29]][,"AWC"] = rowMeans(neu_list2[[29]][,c("AWC_OFF","AWC_ON")])
neu_list2[[44]][,"IL2VL"] = rowMeans(neu_list2[[44]][,c("IL2_DV","IL2_LR")])
neu_list2[[60]][,"RMDDL"] = rowMeans(neu_list2[[60]][,c("RMD_DV","RMD_LR")])
neu_list2[[61]][,"RMDVL"] = rowMeans(neu_list2[[61]][,c("RMD_DV","RMD_LR")])
neu_list2[[29]] = select(neu_list2[[29]],-c("AWC_OFF","AWC_ON"))
neu_list2[[44]] = select(neu_list2[[44]],-c("IL2_DV","IL2_LR"))
neu_list2[[60]] = select(neu_list2[[60]],-c("RMD_DV","RMD_LR"))
neu_list2[[61]] = select(neu_list2[[61]],-c("RMD_DV","RMD_LR"))

#zscore
neu_list3 = list()
for(i in 1:length(neu_list2)){
  neu_list2[[i]] = log2(neu_list2[[i]] + 1)
  neu_list3[[i]] = t(apply(neu_list2[[i]],1,zscore))
}

neu_cell$time_range = paste("X",neu_cell$Start_time,neu_cell$End_time,sep = "_")
for(i in 1:length(neu_list3)){
  for(j in 1:ncol(neu_list3[[i]])){
    for(m in 1:nrow(neu_cell)){
      if(colnames(neu_list3[[i]])[j] == neu_cell[m,"Precursor2"]){
        colnames(neu_list3[[i]])[j] = neu_cell[m,"time_range"]
      }
    }
  }
}
colnames(neu_list3[[29]])[colnames(neu_list3[[29]]) == "AWC" ] = "X_2600_3200"
colnames(neu_list3[[44]])[colnames(neu_list3[[44]]) == "IL2VL" ] = "X_2600_3200"
colnames(neu_list3[[60]])[colnames(neu_list3[[60]]) == "RMDDL" ] = "X_2600_3200"
colnames(neu_list3[[61]])[colnames(neu_list3[[61]]) == "RMDVL" ] = "X_2600_3200"

time_neu = c(seq(55,800,by = 10),1800,3000)
neu_list = list()
for(i in 1:length(neu_list3)){
  df = data.frame(matrix(nrow = nrow(neu_list3[[i]]), ncol = length(time_neu) + 1))
  rownames(df) = rownames(neu_list3[[i]])
  colnames(df) = c(time_neu,"celltype")
  df[,"celltype"] = type_neu[i]
  for(j in 1:ncol(neu_list3[[i]])){
    for(m in 1:(ncol(df) -1)){
      if(isTRUE(as.numeric(colnames(df)[m]) > as.numeric(strsplit(colnames(neu_list3[[i]])[j],"_")[[1]][2]) & as.numeric(colnames(df)[m]) < as.numeric(strsplit(colnames(neu_list3[[i]])[j],"_")[[1]][3]))){
        df[,m] = neu_list3[[i]][,j]
      }
    }
  }
  neu_list[[i]] = df
}


rec_z_cls2$WBGene = rownames(rec_z_cls2)
cls_list = list(filter(rec_z_cls2,cls == 2),filter(rec_z_cls2,cls == 3))
cls_list2 = list()
type_cls = c("cls2","cls3")
for(i in 1:length(neu_list)){
  df = data.frame(matrix(nrow = length(type_cls), ncol = ncol(neu_list[[i]])-1))
  rownames(df) = type_cls
  colnames(df) = colnames(neu_list[[i]])[1:(ncol(neu_list[[i]]) -1)]
  for(j in 1:length(cls_list)){
    cur = neu_list[[i]][rownames(neu_list[[i]]) %in% unique(cls_list[[j]]$WBGene),]
    cur = select(cur,-"celltype")
    df[j,] = apply(cur,2,mean2)
  }
  df$cls_type = rownames(df)
  df$celltype = type_neu[i]
  cls_list2[[i]] = df
}
cls_list2_2 = do.call(rbind,cls_list2)
my_palette <- colorRampPalette(c("lightgreen", "black", "red"))(n = 299)
col_breaks = c(seq(-1.2,-0.3,length=100),  # for green
               seq(-0.31,0.6,length=100),           # for black
               seq(0.61,1.5,length=100))

for(i in 1:length(type_cls)){
  cur = filter(cls_list2_2,cls_type == type_cls[i])
  rownames(cur) = cur$celltype
  cur = select(cur,-c("cls_type","celltype"))
  cur2 = as.matrix(cur)
  
  #get rowder rownames
  aaa3 = heatmap.2(cur2)
  name_ord = rownames(cur2)[aaa3$rowInd]
  cur3 = cur2[rev(name_ord),]
  
  #get lineage name
  cur3 = data.frame(cur3) #note cur3 above is a mtrix
  for(m in 1:nrow(cur3)){
    cur_neu = filter(neu2,Neuron == rownames(cur3)[m])
    for(n in 1:nrow(cur_neu)){
      if(unique(cur_neu$Neuron) == "DVC" & isTRUE(substring(cur_neu[n,"Precursor"],1,1) == "C")){
        cur_neu[n,"num"] = nchar(cur_neu[n,"Precursor"])
      } else if(isTRUE(substring(cur_neu[n,"Precursor"],1,2) == "AB")){
        cur_neu[n,"num"] = nchar(cur_neu[n,"Precursor"])
      }
    }
    cur3[m,"name_lin"] = paste0(rownames(cur3)[m],"@@@",cur_neu[which.max(cur_neu[,"num"]),"Precursor"])
  }
  for(j in 1:ncol(cur3)){
    colnames(cur3)[j] = strsplit(colnames(cur3)[j],"X")[[1]][2]
  }
  
  pdf(paste0("neuron_",type_cls[i],".pdf"),height = 15,width=25)
  heatmap.2(cur2,
            density.info="none",  # turns off density plot inside color legend
            key=T, keysize=1.0,
            scale="none",
            trace="none",         # turns off trace lines inside the heat map
            col=my_palette,       # use on color palette defined earlier
            dendrogram="row",     # only draw a row dendrogram
            Colv="NA",
            srtCol=22,            # turn off column clustering
            margins =c(11,11),
            na.rm = T)
  dev.off()
}

#non-neuronal cells
type_non_neu = sort(unique(non_neu_cell$Lineage))
non_neu_list2 = list()
for(i in 1:length(type_non_neu)){
  for(j in 1:nrow(celltype)){
    cur = filter(non_neu_cell,Lineage == type_non_neu[i])
    cur_type = celltype[celltype$celltype2 %in% cur$Precursor2,]
    cur_exp = select(sc_exp,c(cur_type[,"index"]))
    for(g in 1:ncol(cur_exp)){
      for(h in 1:nrow(cur_type)){
        if(isTRUE(colnames(cur_exp)[g] == cur_type[h,"index"])){
          colnames(cur_exp)[g] = cur_type[h,"celltype2"]
        }}}
  }
  non_neu_list2[[i]] = cur_exp
  if(i %% 10 ==0){
    cat(i,"rows finished \n")
  }
}

#zscore
non_neu_list3 = list()
for(i in 1:length(non_neu_list2)){
  non_neu_list2[[i]] = log2(non_neu_list2[[i]] + 1)
  non_neu_list3[[i]] = t(apply(non_neu_list2[[i]],1,zscore))
}

non_neu_cell$time_range = paste("X",non_neu_cell$Start_time,non_neu_cell$End_time,sep = "_")
for(i in 1:length(non_neu_list3)){
  for(j in 1:ncol(non_neu_list3[[i]])){
    for(m in 1:nrow(non_neu_cell)){
      if(colnames(non_neu_list3[[i]])[j] == non_neu_cell[m,"Precursor2"]){
        colnames(non_neu_list3[[i]])[j] = non_neu_cell[m,"time_range"]
      }
    }
  }
}

time_non_neu = c(seq(55,800,by = 10),1800,3000)
non_neu_list = list()
for(i in 1:length(non_neu_list3)){
  df = data.frame(matrix(nrow = nrow(non_neu_list3[[i]]), ncol = length(time_non_neu) + 1))
  rownames(df) = rownames(non_neu_list3[[i]])
  colnames(df) = c(time_non_neu,"celltype")
  df[,"celltype"] = type_non_neu[i]
  for(j in 1:ncol(non_neu_list3[[i]])){
    for(m in 1:(ncol(df) -1)){
      if(isTRUE(as.numeric(colnames(df)[m]) > as.numeric(strsplit(colnames(non_neu_list3[[i]])[j],"_")[[1]][2]) & as.numeric(colnames(df)[m]) < as.numeric(strsplit(colnames(non_neu_list3[[i]])[j],"_")[[1]][3]))){
        df[,m] = non_neu_list3[[i]][,j]
      }
    }
  }
  non_neu_list[[i]] = df
}

rec_z_cls2$WBGene = rownames(rec_z_cls2)
cls_list = list(filter(rec_z_cls2,cls == 2),filter(rec_z_cls2,cls == 3))
cls_list2 = list()
type_cls = c("cls2","cls3")
for(i in 1:length(non_neu_list)){
  df = data.frame(matrix(nrow = length(type_cls), ncol = ncol(non_neu_list[[i]])-1))
  rownames(df) = type_cls
  colnames(df) = colnames(non_neu_list[[i]])[1:ncol(non_neu_list[[i]]) -1]
  for(j in 1:length(cls_list)){
    cur = non_neu_list[[i]][rownames(non_neu_list[[i]]) %in% unique(cls_list[[j]]$WBGene),]
    cur = select(cur,-"celltype")
    df[j,] = apply(cur,2,mean2)
  }
  df$cls_type = rownames(df)
  df$celltype = type_non_neu[i]
  cls_list2[[i]] = df
}

cls_list2_2 = do.call(rbind,cls_list2)
my_palette <- colorRampPalette(c("lightgreen", "black", "red"))(n = 299)
col_breaks = c(seq(-1.2,-0.3,length=100),  # for green
               seq(-0.31,0.6,length=100),           # for black
               seq(0.61,1.5,length=100))

for(i in 1:length(type_cls)){
  cur = filter(cls_list2_2,cls_type == type_cls[i])
  rownames(cur) = cur$celltype
  cur = select(cur,-c("cls_type","celltype"))
  cur2 = as.matrix(cur)
  
  #get rowder rownames
  aaa3 = heatmap.2(cur2)
  name_ord = rownames(cur2)[aaa3$rowInd]
  cur3 = cur2[rev(name_ord),]
  
  #get lineage name
  cur3 = data.frame(cur3) #note cur3 above is a mtrix
  for(m in 1:nrow(cur3)){
    cur_neu = filter(non_neu_cell,Lineage == rownames(cur3)[m])
    for(n in 1:nrow(cur_neu)){
      cur_neu = cur_neu[cur_neu$Precursor %in% lin_cells$X,]
      cur_neu[n,"num"] = nchar(cur_neu[n,"Precursor"])
    }
    cur3[m,"name_lin"] = paste0(rownames(cur3)[m],"@@@",cur_neu[which.max(cur_neu[,"num"]),"Precursor"])
  }
  for(j in 1:ncol(cur3)){
    colnames(cur3)[j] = strsplit(colnames(cur3)[j],"X")[[1]][2]
  }
  
  pdf(paste0("non_neuron_",type_cls[i],".pdf"),height = 15,width=25)
  heatmap.2(cur2,
            density.info="none",  # turns off density plot inside color legend
            key=T, keysize=1.0,
            scale="none",
            trace="none",         # turns off trace lines inside the heat map
            col=my_palette,       # use on color palette defined earlier
            dendrogram="row",     # only draw a row dendrogram
            Colv="NA",
            srtCol=22,            # turn off column clustering
            margins =c(11,11),
            na.rm = T)
  dev.off()
}
