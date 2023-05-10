library(dplyr)
library(RColorBrewer)
library(gplots)
load("cluster_1_2_3.RData")
load("neuronal_non neuronal_tissues.RData")

#neurons
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
    cur_neu = filter(non_neu_filter,Lineage == rownames(cur3)[m])
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
