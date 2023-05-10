library(dplyr)
library(stringr)
library(ape)
library(hash)

orthogroup <- read.csv("orthogroups.tsv",sep="\t",stringsAsFactors = F)


rapid <- readLines("Expanded_contracted_OGs.txt")
split<- c()
rapidlist <- list()
for(i in 1:length(rapid)){
  ln <-str_match_all(rapid[i],"OG[0-9]+\\[[-~][0-9]+\\*\\]")[[1]]  
  rapidlist[[i]] <- ln[grepl("-",ln,fixed)]
  split <- append(split ,ln[grepl("-",ln,fixed)])
}
for(i in 1:length(split)){
  split[i] <- str_match(split[i] , "OG[0-9]+")[[1]][1]
}

for(i in 1:length(rapidlist)){
  if(length(rapidlist[[i]]) != 0){
    rapidlist[[i]] <- unlist(str_match_all(rapidlist[[i]] , "OG[0-9]+"))
  }else{
    rapidlist[[i]] <- "noorthogroup"
  }
}

expandedOG <- split


#parse result
load("gene_mapping.RData")
outhead <- c("domains" , unique(gene_mapping$species))
pfam_all <- read.csv("pfam_11_species_longest.csv")
resultlist <- list()
genecnt <- data.frame( OGcnt= rep(0,length(split)) , pfamcnt = rep(0,length(split)))

for(i in 1:length(expandedOG)){
  genes <- filter(orthogroup , Orthogroup == expandedOG[i])
  genes <- strsplit( paste(genes[2:12] , collapse = ",") , ",")[[1]]
  genes <- genes[genes!=""]
  genes <- gsub(" +", "",genes)
  curpfam <- filter( pfam_all , name %in% genes)
  genecnt[i,1] <- length(genes)
  genecnt[i,2] <- nrow(curpfam)
  curpfam <- select(curpfam , c(name, domain ,species ))
  if(nrow(curpfam)> 0){
    curpfam$gp <- expandedOG[i]
    colnames(curpfam) <- c("gene" , "domain" , "species","gp")
  }
  
  resultlist[[i]] <- curpfam
}

resultlist2 <- resultlist[lengths(resultlist)!= 3]
domlistname <- expandedOG[lengths(resultlist)!= 3]

orthogroup <- filter(orthogroup , Orthogroup %in% expandedOG)
orthogroup$genecnt <- 0
for(i in 1:nrow(orthogroup)){
  gene <- paste(orthogroup[i,2:12],collapse = ",")
  gene <- length(strsplit(gene , ",")[[1]])
  orthogroup[i,"genecnt"] <- gene
}


for( i in 1:length(resultlist2)){
  x <- resultlist2[[i]]
  x$remove <- F
  x <- distinct(x , gene , domain, species,.keep_all=T)
  domscnt <- table(x$domain)
  names(domscnt) -> doms
  names(domscnt) <- NULL
  gp <- strsplit(x$gp[1],"_")[[1]][1]
  x$gp <- gp
  genecnt <- filter(orthogroup , Orthogroup == gp)$genecnt
  for(j in 1:length(doms)){
   if(domscnt[j]/genecnt < 0.5 ){
     x[x$domain==doms[j] , "remove" ] <- T
   } 
  }
  x <- filter(x , remove ==F) %>% select(-c(remove))
  resultlist2[[i]] <- x
}

#domain name unique
d <- c()
for(i in 1:length(resultlist2)){
  d <- append(d, unique(resultlist2[[i]]$domain))
}
d <- unique(d)

combined_domain <- c("GPCR_combined" , "Fbox_combined" , "BTB_combined" , "HTH_Tnp_Tc3_combined" , "AAA_combined" , "Peptidase_M13_combined", "Histone_combined","Glyco_transf_combined", "Glyco_hydro_combined" , "NHR")

domain_sig_orthogroup <- as.data.frame(matrix(nrow=length(d)+length(combined_domain) , ncol= 2,data=rep("", 2*(length(d)+length(combined_domain))) ))
colnames(domain_sig_orthogroup) <- c("domains","sig_orthogroups")
domain_sig_orthogroup$domains <- c(d , combined_domain )

combined_list <- list()

combined_list[[1]] <- c (d[grepl("7TM|Sre|Srg",d)] ,"srg")
combined_list[[2]] <- c("F-box","FBA_2","FTH","HTH-48")
combined_list[[3]] <- c("BTB","BTB_2")
combined_list[[10]] <- c("zf-C4","Hormone_recep")


result <- as.data.frame(matrix(  ncol=12,nrow=length(d)+length(combined_domain) , data=rep(0 , 12*(length(d)+length(combined_domain)))))
colnames(result)<- outhead
result$domains <- c(d , combined_domain )
sig <- result

for(i in 1:length(resultlist2)){
  if(nrow(resultlist2[[i]])==0){
    next()
  }
  #curgenedom <- as.data.frame(matrix(ncol= 3 , nrow= nrow(resultlist2[[i]]) , data=rep(0,3*nrow(resultlist2[[i]]) )))
  for(j in 1:nrow(resultlist2[[i]])){
      curdomain = resultlist2[[i]][j,"domain"]
      x <- result$domains==curdomain
      y <- which(outhead == resultlist2[[i]][j,"species"])  
      z <- resultlist2[[i]][j,"gp"]
    #check sig. OG per species
      result[ x ,y ] <- paste(result[x,y] , resultlist2[[i]][j,"gene"] ,sep="~")
      domain_sig_orthogroup[  domain_sig_orthogroup$domains ==curdomain ,  2 ] <- paste(domain_sig_orthogroup[  domain_sig_orthogroup$domains == curdomain ,  2 ] ,z ,sep="_")
      
      #check if in combined domain
      for(k in 1:length(combined_list)){
        if(curdomain %in% combined_list[[k]]){
          result[ result$domains==combined_domain[k] , y ] <- paste(result[ result$domains==combined_domain[k],y],resultlist2[[i]][j,"gene"], sep="~")  
        }
      }
    if(z %in% rapidlist[[y-1]]){
      sig[x,y] <- 1
    }
  }
}

for(i in 1:nrow(result)){
  for(j in 2:ncol(result)){
    if(result[i,j]!=0){
    genes <- strsplit( result[i , j ] , "~" ,fixed = T)[[1]]
    genes <- unique(genes[genes!="0"])
    result[ i , j ] <- length(genes)
    }else{
      result[i,j] = 0
    }
    if(sig[i,j]==1){
      result[i,j] <- paste0(result[i,j] , "*")
    }
  }
}

for(i in 1:length(combined_list)){
  for(j in 2:12){
    if(sum(sig[sig$domains %in% combined_list[[i]] , j]) > 0){
      result[result$domains==combined_domain[i] , j ] <- paste0(result[result$domains==combined_domain[i] , j ]  , "*")
    }
  }
}

for(i in 1:length(combined_list)){
  domain_sig_orthogroup[domain_sig_orthogroup$domains==combined_domain[i] , 2 ] <- paste(domain_sig_orthogroup[domain_sig_orthogroup$domains %in% combined_list[[i]] , 2] , collapse = "_")
}

domain_sig_orthogroup$OGcnt <- 0
for(i in 1:nrow(domain_sig_orthogroup)){
  ogs <- unique(strsplit(  domain_sig_orthogroup[i,2] , "_")[[1]])
  ogs <- ogs[ogs!=""]
  domain_sig_orthogroup[i,2] <- paste(ogs , collapse="_")
  domain_sig_orthogroup[i,3] <- length(ogs)
}
domain_sig_orthogroup <- domain_sig_orthogroup[domain_sig_orthogroup$OGcnt!=0,]

# for contracted only
for(i in 1:nrow(result)){
  curog <- strsplit( domain_sig_orthogroup[i,2] , "_")[[1]]
  for(j in 1:length(curog)){
    for(k in 1:11){
      if( curog[j] %in% rapidlist[[k]] ){
        if(result[i , k+1] == "0"){
          result[i , k+1] = "0*"
        }
      }
    } 
  }
}



write.csv(result,"ontracted_domain_total_zerosig.csv", row.names = F)

write.csv(domain_sig_orthogroup , "contracted_domain_OG_number.csv",row.names = F)

```