library(ape)
#doanloading eleven Caenorhabditis species from WormBase Parasite website
species = list.files("path",pattern="fa")#eleven species protein file
fa <- list()

for(i in 1:length(species)){
  fa[[i]]<- as.character(read.FASTA(paste0("path",species[i]),type="AA" ))
}
fanames <- list()
fanames2 <- list()
for(i in 1:17){
  fanames2[[i]] <- fanames[[i]] <- names(fa[[i]])
  for(j in 1:length(fanames[[i]])){
    fanames2[[i]][j] <- strsplit(fanames[[i]][j], " ")[[1]][3]
  }
}
newfa <- list()
for(i in 1:length(fanames)){
  uniquegene <- unique(fanames2[[i]])
  newfa[[i]] <- list()
  newname <- c()
  for( j in 1:length(uniquegene)){
    transcripts <- fa[[i]][ fanames2[[i]] == uniquegene[j]  ] 
    if(length(transcripts)==1){
      newfa[[i]][j] <- transcripts
      newname <- append(newname , names(transcripts))
    }else{
      max = nchar(paste0(transcripts[1],collapse = ""))
      maxk = 1
      for(k in 2:length(transcripts)){
        cur = nchar(paste0(transcripts[k],collapse = ""))
        if(cur > max ){
          max = cur
          maxk = k
        }
      }
      newfa[[i]][j] <- transcripts[maxk]
      newname   <- append(newname , names(transcripts[maxk]))
    }
  }
  names(newfa[[i]]) <- newname
}
species <- gsub(".fa","",species)
for(i in 1:length(species)){
  #newfa[[i]] <- as.AAbin(newfa[[i]])
  write.FASTA(newfa[[i]] , file=paste0("path/longest_transcript/",species[i],"_longest_transcript.fa" )  )
}
