load("duplications.RData")
library(dplyr)
duplication=filter(duplication_separated,Support >= 0.5)
domain=read.csv("Domain_11_species.csv",header=T)
library(dplyr)
#elegans
duplication_elegans=filter(duplication,Species=="elegans")
duplication_N0=duplication_elegans[grep("N0",duplication_elegans$Species.Tree.Node),]
duplication_N3N1=duplication_elegans[grep("N3|N1",duplication_elegans$Species.Tree.Node),]
duplication_N4Cel=duplication_elegans[grep("N4|elegans",duplication_elegans$Species.Tree.Node),]
duplication_N3N1_2=duplication_N3N1[duplication_N3N1$Gene %in% setdiff(duplication_N3N1$Gene, duplication_N4Cel$Gene),]
a=rbind(duplication_N3N1,duplication_N4Cel)
duplication_N0_2=duplication_N0[duplication_N0$Gene %in% setdiff(duplication_N0$Gene,a$Gene),]
domain_elegans=filter(domain,species=="elegans")
domain_elegans_60=data.frame(summary(domain_elegans[domain_elegans$Wb %in% unique(duplication_N4Cel$Gene),]$domain,maxsum=max(lengths(lapply(domain_elegans[domain_elegans$Wb %in% unique(duplication_N4Cel$Gene),],unique)))))
domain_elegans_60_130=data.frame(summary(domain_elegans[domain_elegans$Wb %in% unique(duplication_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_elegans[domain_elegans$Wb %in% unique(duplication_N3N1_2$Gene),],unique)))))
domain_elegans_130=data.frame(summary(domain_elegans[domain_elegans$Wb %in% unique(duplication_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_elegans[domain_elegans$Wb %in% unique(duplication_N0_2$Gene),],unique)))))

#inopinata
duplication_inopinata=filter(duplication,Species=="inopinata")
inopinata_N0=duplication_inopinata[grep("N0",duplication_inopinata$Species.Tree.Node),]
inopinata_N3N1=duplication_inopinata[grep("N3|N1",duplication_inopinata$Species.Tree.Node),]
inopinata_N4Cel=duplication_inopinata[grep("N4|inopinata",duplication_inopinata$Species.Tree.Node),]
inopinata_N3N1_2=inopinata_N3N1[inopinata_N3N1$Gene %in% setdiff(inopinata_N3N1$Gene, inopinata_N4Cel$Gene),]
a=rbind(inopinata_N3N1,inopinata_N4Cel)
inopinata_N0_2=inopinata_N0[inopinata_N0$Gene %in% setdiff(inopinata_N0$Gene,a$Gene),]
domain_inopinata=filter(domain,species=="inopinata")
domain_inopinata_60=data.frame(summary(domain_inopinata[domain_inopinata$name %in% unique(inopinata_N4Cel$Gene),]$domain,maxsum=max(lengths(lapply(domain_inopinata[domain_inopinata$name %in% unique(inopinata_N4Cel$Gene),],unique)))))
domain_inopinata_60_130=data.frame(summary(domain_inopinata[domain_inopinata$name %in% unique(inopinata_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_inopinata[domain_inopinata$name %in% unique(inopinata_N3N1_2$Gene),],unique)))))
domain_inopinata_130=data.frame(summary(domain_inopinata[domain_inopinata$name %in% unique(inopinata_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_inopinata[domain_inopinata$name %in% unique(inopinata_N0_2$Gene),],unique)))))

#bovis
duplication_bovis=filter(duplication,Species=="bovis")
bovis_N0=duplication_bovis[grep("N0",duplication_bovis$Species.Tree.Node),]
bovis_bov=duplication_bovis[grep("bovis",duplication_bovis$Species.Tree.Node),]
bovis_N0_2=bovis_N0[bovis_N0$Gene %in% setdiff(bovis_N0$Gene, bovis_bov$Gene),]
domain_bovis=filter(domain,species=="bovis")
domain_bovis_in130=data.frame(summary(domain_bovis[domain_bovis$name %in% unique(bovis_bov$Gene),]$domain,maxsum=max(lengths(lapply(domain_bovis[domain_bovis$name %in% unique(bovis_bov$Gene),],unique)))))
domain_bovis_before130=data.frame(summary(domain_bovis[domain_bovis$name %in% unique(bovis_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_bovis[domain_bovis$name %in% unique(bovis_N0_2$Gene),],unique)))))

#becei
duplication_becei=filter(duplication,Species=="becei")
becei_N0=duplication_becei[grep("N0",duplication_becei$Species.Tree.Node),]
becei_N1=duplication_becei[grep("N1",duplication_becei$Species.Tree.Node),]
becei_N2bec=duplication_becei[grep("N2|becei",duplication_becei$Species.Tree.Node),]
becei_N1_2=becei_N1[becei_N1$Gene %in% setdiff(becei_N1$Gene, becei_N2bec$Gene),]
a=rbind(becei_N1,becei_N2bec)
becei_N0_2=becei_N0[becei_N0$Gene %in% setdiff(becei_N0$Gene,a$Gene),]
domain_becei=filter(domain,species=="becei")
domain_becei_60=data.frame(summary(domain_becei[domain_becei$name %in% unique(becei_N2bec$Gene),]$domain,maxsum=max(lengths(lapply(domain_becei[domain_becei$name %in% unique(becei_N2bec$Gene),],unique)))))
domain_becei_60_130=data.frame(summary(domain_becei[domain_becei$name %in% unique(becei_N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_becei[domain_becei$name %in% unique(becei_N1_2$Gene),],unique)))))
domain_becei_130=data.frame(summary(domain_becei[domain_becei$name %in% unique(becei_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_becei[domain_becei$name %in% unique(becei_N0_2$Gene),],unique)))))

#panamensis
duplication_panamensis=filter(duplication,Species=="panamensis")
panamensis_N0=duplication_panamensis[grep("N0",duplication_panamensis$Species.Tree.Node),]
panamensis_N1=duplication_panamensis[grep("N1",duplication_panamensis$Species.Tree.Node),]
panamensis_N2pan=duplication_panamensis[grep("N2|panamensis",duplication_panamensis$Species.Tree.Node),]
panamensis_N1_2=panamensis_N1[panamensis_N1$Gene %in% setdiff(panamensis_N1$Gene, panamensis_N2pan$Gene),]
a=rbind(panamensis_N1,panamensis_N2pan)
panamensis_N0_2=panamensis_N0[panamensis_N0$Gene %in% setdiff(panamensis_N0$Gene,a$Gene),]
domain_panamensis=filter(domain,species=="panamensis")
domain_panamensis_60=data.frame(summary(domain_panamensis[domain_panamensis$name %in% unique(panamensis_N2pan$Gene),]$domain,maxsum=max(lengths(lapply(domain_panamensis[domain_panamensis$name %in% unique(panamensis_N2pan$Gene),],unique)))))
domain_panamensis_60_130=data.frame(summary(domain_panamensis[domain_panamensis$name %in% unique(panamensis_N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_panamensis[domain_panamensis$name %in% unique(panamensis_N1_2$Gene),],unique)))))
domain_panamensis_130=data.frame(summary(domain_panamensis[domain_panamensis$name %in% unique(panamensis_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_panamensis[domain_panamensis$name %in% unique(panamensis_N0_2$Gene),],unique)))))

#tropicalis
duplication_tropicalis=filter(duplication,Species=="tropicalis")
tropicalis_N0=duplication_tropicalis[grep("N0",duplication_tropicalis$Species.Tree.Node),]
tropicalis_N3N1=duplication_tropicalis[grep("N3|N1",duplication_tropicalis$Species.Tree.Node),]
tropicalis_N5tro=duplication_tropicalis[grep("N5|tropicalis",duplication_tropicalis$Species.Tree.Node),]
tropicalis_N3N1_2=tropicalis_N3N1[tropicalis_N3N1$Gene %in% setdiff(tropicalis_N3N1$Gene, tropicalis_N5tro$Gene),]
a=rbind(tropicalis_N3N1,tropicalis_N5tro)
tropicalis_N0_2=tropicalis_N0[tropicalis_N0$Gene %in% setdiff(tropicalis_N0$Gene,a$Gene),]
domain_tropicalis=filter(domain,species=="tropicalis")
domain_tropicalis_60=data.frame(summary(domain_tropicalis[domain_tropicalis$name %in% unique(tropicalis_N5tro$Gene),]$domain,maxsum=max(lengths(lapply(domain_tropicalis[domain_tropicalis$name %in% unique(tropicalis_N5tro$Gene),],unique)))))
domain_tropicalis_60_130=data.frame(summary(domain_tropicalis[domain_tropicalis$name %in% unique(tropicalis_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_tropicalis[domain_tropicalis$name %in% unique(tropicalis_N3N1_2$Gene),],unique)))))
domain_tropicalis_130=data.frame(summary(domain_tropicalis[domain_tropicalis$name %in% unique(tropicalis_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_tropicalis[domain_tropicalis$name %in% unique(tropicalis_N0_2$Gene),],unique)))))

#remanei
duplication_remanei=filter(duplication,Species=="remanei")
remanei_N0=duplication_remanei[grep("N0",duplication_remanei$Species.Tree.Node),]
remanei_N3N1=duplication_remanei[grep("N3|N1",duplication_remanei$Species.Tree.Node),]
remanei_N8N6N5rem=duplication_remanei[grep("N8|N6|N5|remanei",duplication_remanei$Species.Tree.Node),]
remanei_N3N1_2=remanei_N3N1[remanei_N3N1$Gene %in% setdiff(remanei_N3N1$Gene, remanei_N8N6N5rem$Gene),]
a=rbind(remanei_N3N1,remanei_N8N6N5rem)
remanei_N0_2=remanei_N0[remanei_N0$Gene %in% setdiff(remanei_N0$Gene,a$Gene),]
domain_remanei=filter(domain,species=="remanei")
domain_remanei_60=data.frame(summary(domain_remanei[domain_remanei$name %in% unique(remanei_N8N6N5rem$Gene),]$domain,maxsum=max(lengths(lapply(domain_remanei[domain_remanei$name %in% unique(remanei_N8N6N5rem$Gene),],unique)))))
domain_remanei_60_130=data.frame(summary(domain_remanei[domain_remanei$name %in% unique(remanei_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_remanei[domain_remanei$name %in% unique(remanei_N3N1_2$Gene),],unique)))))
domain_remanei_130=data.frame(summary(domain_remanei[domain_remanei$name %in% unique(remanei_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_remanei[domain_remanei$name %in% unique(remanei_N0_2$Gene),],unique)))))

#latens
duplication_latens=filter(duplication,Species=="latens")
latens_N0=duplication_latens[grep("N0",duplication_latens$Species.Tree.Node),]
latens_N3N1=duplication_latens[grep("N3|N1",duplication_latens$Species.Tree.Node),]
latens_N8N6N5lat=duplication_latens[grep("N8|N6|N5|latens",duplication_latens$Species.Tree.Node),]
latens_N3N1_2=latens_N3N1[latens_N3N1$Gene %in% setdiff(latens_N3N1$Gene, latens_N8N6N5lat$Gene),]
a=rbind(latens_N3N1,latens_N8N6N5lat)
latens_N0_2=latens_N0[latens_N0$Gene %in% setdiff(latens_N0$Gene,a$Gene),]
domain_latens=filter(domain,species=="latens")
domain_latens_60=data.frame(summary(domain_latens[domain_latens$name %in% unique(latens_N8N6N5lat$Gene),]$domain,maxsum=max(lengths(lapply(domain_latens[domain_latens$name %in% unique(latens_N8N6N5lat$Gene),],unique)))))
domain_latens_60_130=data.frame(summary(domain_latens[domain_latens$name %in% unique(latens_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_latens[domain_latens$name %in% unique(latens_N3N1_2$Gene),],unique)))))
domain_latens_130=data.frame(summary(domain_latens[domain_latens$name %in% unique(latens_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_latens[domain_latens$name %in% unique(latens_N0_2$Gene),],unique)))))

#tribulationis
duplication_tribulationis=filter(duplication,Species=="tribulationis")
tribulationis_N0=duplication_tribulationis[grep("N0",duplication_tribulationis$Species.Tree.Node),]
tribulationis_N3N1=duplication_tribulationis[grep("N3|N1",duplication_tribulationis$Species.Tree.Node),]
tribulationis_N7N6N5tri=duplication_tribulationis[grep("N7|N6|N5|tribulationis",duplication_tribulationis$Species.Tree.Node),]
tribulationis_N3N1_2=tribulationis_N3N1[tribulationis_N3N1$Gene %in% setdiff(tribulationis_N3N1$Gene, tribulationis_N7N6N5tri$Gene),]
a=rbind(tribulationis_N3N1,tribulationis_N7N6N5tri)
tribulationis_N0_2=tribulationis_N0[tribulationis_N0$Gene %in% setdiff(tribulationis_N0$Gene,a$Gene),]
domain_tribulationis=filter(domain,species=="tribulationis")
domain_tribulationis_60=data.frame(summary(domain_tribulationis[domain_tribulationis$name %in% unique(tribulationis_N7N6N5tri$Gene),]$domain,maxsum=max(lengths(lapply(domain_tribulationis[domain_tribulationis$name %in% unique(tribulationis_N7N6N5tri$Gene),],unique)))))
domain_tribulationis_60_130=data.frame(summary(domain_tribulationis[domain_tribulationis$name %in% unique(tribulationis_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_tribulationis[domain_tribulationis$name %in% unique(tribulationis_N3N1_2$Gene),],unique)))))
domain_tribulationis_130=data.frame(summary(domain_tribulationis[domain_tribulationis$name %in% unique(tribulationis_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_tribulationis[domain_tribulationis$name %in% unique(tribulationis_N0_2$Gene),],unique)))))

#briggsae
duplication_briggsae=filter(duplication,Species=="briggsae")
briggsae_N0=duplication_briggsae[grep("N0",duplication_briggsae$Species.Tree.Node),]
briggsae_N3N1=duplication_briggsae[grep("N3|N1",duplication_briggsae$Species.Tree.Node),]
briggsae_N9N7N6N5bri=duplication_briggsae[grep("N9|N7|N6|N5|briggsae",duplication_briggsae$Species.Tree.Node),]
briggsae_N3N1_2=briggsae_N3N1[briggsae_N3N1$Gene %in% setdiff(briggsae_N3N1$Gene, briggsae_N9N7N6N5bri$Gene),]
a=rbind(briggsae_N3N1,briggsae_N9N7N6N5bri)
briggsae_N0_2=briggsae_N0[briggsae_N0$Gene %in% setdiff(briggsae_N0$Gene,a$Gene),]
domain_briggsae=filter(domain,species=="briggsae")
domain_briggsae2=merge(domain_briggsae,bri_WB,by="name")
domain_briggsae2_60=data.frame(summary(domain_briggsae2[domain_briggsae2$Gene %in% unique(briggsae_N9N7N6N5bri$Gene),]$domain,maxsum=max(lengths(lapply(domain_briggsae2[domain_briggsae2$Gene %in% unique(briggsae_N9N7N6N5bri$Gene),],unique)))))
domain_briggsae2_60_130=data.frame(summary(domain_briggsae2[domain_briggsae2$Gene %in% unique(briggsae_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_briggsae2[domain_briggsae2$Gene %in% unique(briggsae_N3N1_2$Gene),],unique)))))
domain_briggsae2_130=data.frame(summary(domain_briggsae2[domain_briggsae2$Gene %in% unique(briggsae_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_briggsae2[domain_briggsae2$Gene %in% unique(briggsae_N0_2$Gene),],unique)))))

#nigoni
duplication_nigoni=filter(duplication,Species=="nigoni")
nigoni_N0=duplication_nigoni[grep("N0",duplication_nigoni$Species.Tree.Node),]
nigoni_N3N1=duplication_nigoni[grep("N3|N1",duplication_nigoni$Species.Tree.Node),]
nigoni_N9N7N6N5nig=duplication_nigoni[grep("N9|N7|N6|N5|nigoni",duplication_nigoni$Species.Tree.Node),]
nigoni_N3N1_2=nigoni_N3N1[nigoni_N3N1$Gene %in% setdiff(nigoni_N3N1$Gene, nigoni_N9N7N6N5nig$Gene),]
a=rbind(nigoni_N3N1,nigoni_N9N7N6N5nig)
nigoni_N0_2=nigoni_N0[nigoni_N0$Gene %in% setdiff(nigoni_N0$Gene,a$Gene),]
domain_nigoni=filter(domain,species=="nigoni")
domain_nigoni_60=data.frame(summary(domain_nigoni[domain_nigoni$name %in% unique(nigoni_N9N7N6N5nig$Gene),]$domain,maxsum=max(lengths(lapply(domain_nigoni[domain_nigoni$name %in% unique(nigoni_N9N7N6N5nig$Gene),],unique)))))
domain_nigoni_60_130=data.frame(summary(domain_nigoni[domain_nigoni$name %in% unique(nigoni_N3N1_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_nigoni[domain_nigoni$name %in% unique(nigoni_N3N1_2$Gene),],unique)))))
domain_nigoni_130=data.frame(summary(domain_nigoni[domain_nigoni$name %in% unique(nigoni_N0_2$Gene),]$domain,maxsum=max(lengths(lapply(domain_nigoni[domain_nigoni$name %in% unique(nigoni_N0_2$Gene),],unique)))))


