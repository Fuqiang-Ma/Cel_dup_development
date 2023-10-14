library(stats)
library(dplyr)
library(tidyr)
load("cluster_1_2_3.RData")
load("Cel_dup_N4Cel.RData")
load("Single_copy_genes.RData")
load("bulk_dev_RNAseq.RData")
load("celltype.RData")
load("L4_exp.RData")

rec_z_cls2$WBGene = rownames(rec_z_cls2)
cls1 = filter(rec_z_cls2, cls == 1)
cls2 = filter(rec_z_cls2, cls == 2)
cls3 = filter(rec_z_cls2, cls == 3)

male_female2$fc = as.numeric(male_female2$male)/as.numeric(male_female2$female)
male_female2_fc3 = filter(male_female2,fc>=3)
n4cel_single = list(Cel_dup_N4Cel,genes_1to1,cls1,cls2,cls3)
type= c("N4cel","single_copy","cls1","cls2","cls3")
enr_type = data.frame(type)
for(i in 1:length(type)){
  cur = n4cel_single[[i]][n4cel_single[[i]]$WBGene %in% unique(male_female2_fc3$Sample),]
  enr_type[i,"gene_fc3"] = length(unique(cur$WBGene))
  enr_type[i,"total"] = length(unique(n4cel_single[[i]]$WBGene))
  enr_type[i,"enrich"] = phyper(c(enr_type[i,"gene_fc3"]-1),enr_type[i,"total"],c(19997 - enr_type[i,"total"]),length(unique(male_female2_fc3$Sample)),lower.tail = F)
  enr_type[i,"enrich2"] = -log10(enr_type[i,"enrich"])
}

#correlation
orthogroup <- read.csv("./Orthogroups.tsv",sep="\t")
orthogroup3 = select(orthogroup,c("Orthogroup","elegans_longest_transcript"))
names(orthogroup3)[names(orthogroup3) == "elegans_longest_transcript"] = "WBGene"
gene_count = read.csv("Orthogroups.GeneCount.tsv",sep = "\t")
gene_count = select(gene_count,c("Orthogroup","elegans_longest_transcript"))
orthogroup3 = merge(orthogroup3,gene_count,by="Orthogroup")
for(i in 1:nrow(bulk_RNAseq)){
  bulk_RNAseq[i,"max"] = max(as.numeric(bulk_RNAseq[i,2:49]))
}
bulk_RNAseq2 = filter(bulk_RNAseq,max>10)
fpkm_less10 = Genes19997[Genes19997$WBGene %in% setdiff(Genes19997$WBGene,bulk_RNAseq2$WBGene),]
rec_z_cls2$WBGene = rownames(rec_z_cls2)
fpkm_others = bulk_RNAseq2[bulk_RNAseq2$WBGene %in% setdiff(bulk_RNAseq2$WBGene,rec_z_cls2$WBGene),]
cls_fpkm = list(cls1,cls2,cls3,fpkm_less10,fpkm_others)
type6 = c("cls1","cls2","cls3","fpkm_less10","fpkm_others")
for(i in 1:nrow(orthogroup3)){
  if(orthogroup3[i,"WBGene"] != ""){
    cur = strsplit(orthogroup3[i,"WBGene"],",")[[1]]
    for(j in 1:length(cls_fpkm)){
      orthogroup3[i,type6[j]] = length(intersect(cls_fpkm[[j]]$WBGene,cur))
      orthogroup3[i,paste0(type6[j],"_perc")] = length(intersect(cls_fpkm[[j]]$WBGene,cur))/length(cur)
    }
  }
  if(i %% 8000 == 0){
    cat(i,"rows finished \n")
  }
}
orthogroup3_2 = filter(orthogroup3, cls1 >0 | cls2 >0 | cls3 >0)
sum(filter(orthogroup3,cls1 > 1 & cls2 ==0 & cls3 == 0)$cls1)
sum(filter(orthogroup3,cls2 > 1 & cls1 ==0 & cls3 == 0)$cls2)
sum(filter(orthogroup3,cls3 > 1 & cls2 ==0 & cls1 == 0)$cls3)

zscore = function(x){
  x = (x - mean(x))/sd(x)
}

aaa2 = celltype[celltype$celltype2 %in% colnames(L4_exp),]
L4_tpm = sc_exp[,aaa2$index]
for(i in 1:ncol(L4_tpm)){
  for(j in 1:nrow(aaa2)){
    if(colnames(L4_tpm)[i] == aaa2[j,"index"]){
      colnames(L4_tpm)[i] = aaa2[j,"celltype2"]
    }
  }}
L4_zscore = data.frame(t(apply(log2(L4_tpm + 1),1,zscore)))
L4_tpm$WBGene = rownames(L4_tpm)
cls3_ortho = filter(orthogroup3_2,cls3 >1)
cls3_ortho2 = select(cls3_ortho,c("Orthogroup","WBGene","cls3","cls3_perc"))
cls3_ortho2 = separate_rows(cls3_ortho2,WBGene)
cls3_ortho2 = cls3_ortho2[cls3_ortho2$WBGene %in% intersect(cls3$WBGene,cls3_ortho2$WBGene),]
cls3_ortho_L4 = merge(cls3_ortho2,L4_tpm,by= "WBGene")

bbb = unique(cls3_ortho$Orthogroup)
cur_L4 = merge(cls3_ortho2,L4_tpm,by="WBGene")
df=list()
for (j in 1:length(bbb)){
  df[[j]]=cur_L4[cur_L4$Orthogroup %in% bbb[j],]
}
df[which(lapply(df,nrow) >= 2)] -> x
bbb[which(lapply(df,nrow) >= 2)] -> bbb
for (k in 1:length(bbb)){
  rownames(x[[k]])=x[[k]]$WBGene
  x[[k]]=x[[k]][,c(5:ncol(x[[k]]))]}
expand20=list()
expand20p=list()
for (m in 1:length(bbb)){
  expand20[[m]]=c()
  expand20p[[m]]=c()
  cur = as.data.frame(x[[m]])
  cnt=nrow(cur)
  index=1
  for (h in 1:(cnt-1)){
    for (o in (h+1):cnt){
      if (sum(cur[o,]) == 0){
        next()}
      curtest=cor.test(unlist(cur[h,]),unlist(cur[o,]))
      curcor= as.double(curtest$estimate)
      curcorp= as.double(curtest$p.value)
      if(index ==1){
        expand20[[m]]=curcor
        expand20p[[m]]=curcorp
        index=2} else{
          expand20[[m]]=append(expand20[[m]],curcor)
          expand20p[[m]]=append(expand20p[[m]],curcorp)}
    }}}

corhist_cls3 <- as.data.frame(matrix(ncol=3))
for(n in 1:length(bbb)){
  if (length(expand20[[n]]) == 0){
    next()}
  cur <- as.data.frame(matrix(ncol=3,nrow=length(expand20[[n]])))
  cur[,2] <- expand20[[n]]
  cur[,3] <- expand20p[[n]]
  cur[,1] <- nrow(x[[n]])
  cur[,4] <- bbb[n]
  corhist_cls3 <- bind_rows(corhist_cls3, cur)
}
colnames(corhist_cls3) <- c("Gene_number_in_this_OG","cor" , "pvalue","Orthogroup")
corhist_cls3$Type = "cls2"
corhist_cls3 <- corhist_cls3[2:nrow(corhist_cls3),]

#single copy genes
genes1to1_neuron=merge(genes_1to1,L4_tpm,by = "WBGene")
names(genes1to1_neuron)[names(genes1to1_neuron) == "elegans_longest_transcript"] = "WBGene"
for (i in 1:nrow(genes1to1_neuron)){
  genes1to1_neuron[i,"percent50"]=sum(genes1to1_neuron[i,2:ncol(genes1to1_neuron)]>0,na.rm=T)
}
genes1to1_neuron_less50=filter(genes1to1_neuron,percent50 < 82)

aaa=replicate(915,{x=sample(genes1to1_neuron_less50$WBGene,size=2,replace=F)})
bbb=list()
for (i in 1:915) {
  bbb[[i]]=genes1to1_neuron_less50[genes1to1_neuron_less50$WBGene %in% aaa[,i],]
}

single=list()
singlep=list()
for (i in 1:915){
  index=1
  curtest=cor.test(unlist(bbb[[i]][1,3:c(ncol(bbb[[i]]) -1)]),unlist(bbb[[i]][2,3:c(ncol(bbb[[i]]) -1)]))
  curcor= as.double(curtest$estimate)
  curcorp= as.double(curtest$p.value)
  if(index ==1){
    single[[i]]=curcor
    singlep[[i]]=curcorp
    index=2} else{
      single[[i]]=append(single[[i]],curcor)
      singlep[[i]]=append(singlep[[i]],curcorp)}
}
corhist_single <- as.data.frame(matrix(ncol=3))
for(i in 1:915){
  if (length(single[[i]]) == 0){
    next()}
  cur <- as.data.frame(matrix(ncol=3,nrow=length(single[[i]])))
  cur[,2] <- single[[i]]
  cur[,3] <- singlep[[i]]
  cur[,1] <- nrow(bbb[[i]])
  corhist_single <- bind_rows(corhist_single, cur)
}
colnames(corhist_single) <- c("Gene_number_in_this_OG","cor","pvalue")
corhist_single$Type <- "singleton"
corhist_single <- corhist_single[2:nrow(corhist_single),]
corhist2=bind_rows(corhist_cls3,corhist_single)
pdf("cls3_single_copy.pdf")
plot(ecdf(corhist_cls3$cor),verticals = TRUE,do.points =F,ylab="Cumulative Fraction",xlab="pairwise correlation coefficient of expression",xlim=c(0,1),main="",col="orange",lwd=6,cex.lab=1.5,cex.axis=2)
lines(ecdf(corhist_single$cor),verticals = TRUE, do.points =F,col="gray",lwd=6)
legend("bottomright",legend=c("cor_cls3","cor_single"),lty=c(1,1),col=c("orange","gray"),ncol=1,bty="n",inset=0.03,cex=1.6,lwd=6,text.font = 2)
dev.off()




