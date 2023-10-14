library(stats)
library(purrr)
library(dplyr)
load("bulk_dev_RNAseq.RData")
load("Cel_dup_N4Cel.RData")
load("Single_copy_genes.RData")
load("genes_dNdS.RData")
load("genes_Cel_Cinop_dNdS.RData")
load("cluster_1_2_3.RData")
load("PSG.RData")
load("pi_H_pNpS.RData")
load("divergent_genes.RData")
load("Genes19997.RData")

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
genes_dNdS=filter(genes_dNdS,ds > 0.0005 & dn > 0.0005 & w < 999 & ds <= 2) # excluding genes with ds > 2
genes_dNdS=filter(genes_Cel_Cinop_dNdS,ds > 0.0005 & dn > 0.0005 & w < 999) # for dN/dS by orthologs between C.ele and C.inop

pdf("dNdS_n4cel_others.pdf")
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

###ks test
#n4cel
g1 = list(genes_dNdS[genes_dNdS$wb %in% unique(n4cel_null$WBGene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(n4cel_low$WBGene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(n4cel_medium$WBGene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(n4cel_high$WBGene),]$w)
type = c("n4cel_null","n4cel_low","n4cel_medium","n4cel_high")
alt = combn(1:4,2)
g1_p = data.frame(matrix(nrow = 6, ncol = 3))
colnames(g1_p) = c("type1","type2","pvalue")
for(i in 1:nrow(g1_p)){
  cur1 = alt[,i][1]
  cur2 = alt[,i][2]
  cur = ks.test(g1[[cur1]], g1[[cur2]])
  g1_p[i,"type1"] = type[cur1]
  g1_p[i,"type2"] = type[cur2]
  g1_p[i,"pvalue"] = cur$p.value*6
}

#other genes
g1 = list(genes_dNdS[genes_dNdS$wb %in% unique(genes_null2$WBGene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(genes_low2$WBGene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(genes_medium2$WBGene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(genes_high2$WBGene),]$w)
type = c("others_null","others_low","others_medium","others_high")
alt = combn(1:4,2)
g1_p = data.frame(matrix(nrow = 6, ncol = 3))
colnames(g1_p) = c("type1","type2","pvalue")
for(i in 1:nrow(g1_p)){
  cur1 = alt[,i][1]
  cur2 = alt[,i][2]
  cur = ks.test(g1[[cur1]], g1[[cur2]])
  g1_p[i,"type1"] = type[cur1]
  g1_p[i,"type2"] = type[cur2]
  g1_p[i,"pvalue"] = cur$p.value*6
}

#cluster1, cluster2, cluster3
g1 = list(genes_dNdS$w,genes_dNdS[genes_dNdS$wb %in% rownames(cls1),]$w,
          genes_dNdS[genes_dNdS$wb %in% rownames(cls2),]$w,genes_dNdS[genes_dNdS$wb %in% rownames(cls3),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(duplication_N4Cel399$Gene),]$w,
          genes_dNdS[genes_dNdS$wb %in% unique(genes_1to1$WBGene),]$w)
type = c("all","cls1","cls2","cls3","n4cel","single_copy")
alt = combn(1:6,2)
g1_p = data.frame(matrix(nrow = 15, ncol = 3))
colnames(g1_p) = c("type1","type2","pvalue")
for(i in 1:nrow(g1_p)){
  cur1 = alt[,i][1]
  cur2 = alt[,i][2]
  cur = ks.test(g1[[cur1]], g1[[cur2]])
  g1_p[i,"type1"] = type[cur1]
  g1_p[i,"type2"] = type[cur2]
  g1_p[i,"pvalue"] = cur$p.value*15
}

###hyperdivergent genes
length(unique(divergent_genes2$Wb))
length(intersect(unique(divergent_genes2$Wb),unique(cls1$WBGene)))
length(intersect(unique(divergent_genes2$Wb),unique(cls2$WBGene)))
length(intersect(unique(divergent_genes2$Wb),unique(cls3$WBGene)))

###population genetics statistics excluding hyperdivergent genes
all=list(pi_H_pNpS,genes_1to1,Cel_dup_N4Cel,cls1,cls2,cls3)
genes19997_nohyper = Genes19997[Genes19997$WBGene %in% setdiff(Genes19997$WBGene,unique(divergent_genes2$Wb)),]
all2 = list()
for(i in 1:length(all)){
  all2[[i]] = all[[i]][all[[i]]$WBGene %in% intersect(all[[i]]$WBGene,genes19997_nohyper$WBGene),]
}

#pi
result = as.data.frame(matrix(nrow=3, ncol=6,data=0))
for (i in 1:6){
  a=filter(pi_H_pNpS,Pi >= 0.001)
  result[1,i] = nrow(a[a$WBGene %in% all2[[i]]$WBGene,])
  b=filter(pi_H_pNpS,Pi < 0.001 & Pi >= 0.0001)
  result[2,i] = nrow(b[b$WBGene %in% all2[[i]]$WBGene,])
  c=filter(pi_H_pNpS,Pi < 0.0001)
  result[3,i] = nrow(c[c$WBGene %in% all2[[i]]$WBGene,])
}
colnames(result)=c("All_genes","single_copy_genes","N4Cel","cls1","cls2","cls3")
rownames(result)=c("pi_0.001","pi_0.0001_0.001","pi_0.0001")
result2=as.data.frame(matrix(nrow=3,ncol=6))
for (i in 1:3) {
  for (j in 1:6){
    result2[i,j]=result[i,j]/sum(result[,j])
  }}
colnames(result2)=c("All_genes","single_copy_genes","N4Cel","cls1","cls2","cls3")
rownames(result2)=c("pi_0.001","pi_0.0001_0.001","pi_0.0001")
t(result2)

#pNpS
aa = filter(pi_H_pNpS,pNpS < 10000)
result = as.data.frame(matrix(nrow=3, ncol=6,data=0))
for (i in 1:6){
  a=filter(aa,pNpS >= 1)
  result[1,i] = nrow(a[a$WBGene %in% all2[[i]]$WBGene,])
  b=filter(aa,pNpS < 1 & pNpS >= 0.5)
  result[2,i] = nrow(b[b$WBGene %in% all2[[i]]$WBGene,])
  c=filter(aa,pNpS < 0.5)
  result[3,i] = nrow(c[c$WBGene %in% all2[[i]]$WBGene,])
}
colnames(result)=c("All_genes","single_copy_genes","N4Cel","cls1","cls2","cls3")
rownames(result)=c("pnps_1","pnps_0.5_1","pnps_0.5")
result2=as.data.frame(matrix(nrow=3,ncol=6))
for (i in 1:3) {
  for (j in 1:6){
    result2[i,j]=result[i,j]/sum(result[,j])
  }}
colnames(result2)=c("All_genes","single_copy_genes","N4Cel","cls1","cls2","cls3")
rownames(result2)=c("pnps_1","pnps_0.5_1","pnps_0.5")
t(result2)

#FayWuH
result = as.data.frame(matrix(nrow=4, ncol=6,data=0))
for (i in 1:6){
  a=filter(pi_H_pNpS,FayWuH >= 1)
  result[1,i] = nrow(a[a$WBGene %in% all2[[i]]$WBGene,])
  b=filter(pi_H_pNpS,FayWuH < 1 & FayWuH >= -1)
  result[2,i] = nrow(b[b$WBGene %in% all2[[i]]$WBGene,])
  c=filter(pi_H_pNpS,FayWuH < -1 & FayWuH >= -3)
  result[3,i] = nrow(c[c$WBGene %in% all2[[i]]$WBGene,])
  d=filter(pi_H_pNpS,FayWuH < -3)
  result[4,i] = nrow(d[d$WBGene %in% all2[[i]]$WBGene,])
}
colnames(result)=c("All_genes","single_copy_genes","N4Cel","cls1","cls2","cls3")
rownames(result)=c("FayWuH_1","FayWuH_-1_1","FayWuH_-3_-1","FayWuH_-3")
result2=as.data.frame(matrix(nrow=4,ncol=6))
for (i in 1:4) {
  for (j in 1:6){
    result2[i,j]=result[i,j]/sum(result[,j])
  }}
colnames(result2)=c("All_genes","single_copy_genes","N4Cel","cls1","cls2","cls3")
rownames(result2)=c("FayWuH_1","FayWuH_-1_1","FayWuH_-3_-1","FayWuH_-3")
t(result2)








