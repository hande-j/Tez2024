## MAF Calculation of Imputed vs Called Data April24
library(dplyr)
library(ggplot2)

setwd("/mnt/NEOGENE4/projects/medical_2020/glimpse/sixtyfour/freqs/")

load("/mnt/NEOGENE4/projects/medical_2020/glimpse/sixtyfour/freqs/snpmafcalculation.RData")

kgpmaf5hosnpcount <- dim(read.table(gzfile("/mnt/NEOGENE4/projects/medical_2020/glimpse/fiftysix/KGP_panel/1KGpanel.maf5.ho.sites.tsv.gz")))[1] 

kgpmaf1_5hosnpcount <- dim(read.table(gzfile("/mnt/NEOGENE4/projects/medical_2020/glimpse/fiftysix/KGP_panel/1KGpanel.maf1_5.ho.sites.tsv.gz")))[1]  


snplistmaf5 <- read.table(gzfile("/mnt/NEOGENE4/projects/medical_2020/glimpse/fiftysix/KGP_panel/1KGpanel.maf5.ho.sites.tsv.gz")) #%>% View()

snplistmaf1_5 <- read.table(gzfile("/mnt/NEOGENE4/projects/medical_2020/glimpse/fiftysix/KGP_panel/1KGpanel.maf1_5.ho.sites.tsv.gz")) #%>% View()


set.seed(13)

maf5indexes <- sample(kgpmaf5hosnpcount, 100)

maf1_5indexes <- sample(kgpmaf1_5hosnpcount, 100)


filename=paste("kgpmaf5ho_snplist.txt")
write.table(maf5indexes, file = filename, quote = F, row.names = F, col.names = F )

filename=paste("kgpmaf1_5ho_snplist.txt")
write.table(maf1_5indexes, file = filename, quote = F, row.names = F, col.names = F )

