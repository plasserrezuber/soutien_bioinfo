library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

#localisation donnees entree
dirin<-"Y:/ANALYSES/bin/soutien_bioinfo/BOUGUENNEC/results/"
setwd(dirin)

classi=read_tsv("Y:/ANALYSES/results/wheatomics_wp1/annotation/TE/TE_classification.txt", col_names=T)

d=read_tsv("Length_TE_Lo7_SUPERFamLevel.tsv", col_names=F)
names(d)<-c("chrom", "fam", "length")
ls.str(d)

d=d%>%mutate(genome=str_extract(chrom, "[:upper:]"), chrom=gsub("chr", "", chrom))
d=left_join(d, classi, by="fam")


total_dna=6206789216
d_genome=d%>%filter(chrom!="Un")%>%group_by(fam,fullname)%>%summarise(cumul_length=sum(length), percent_DNA=round(sum(length)/total_dna*100,2))
d_genome=d_genome%>%cbind(genome="R")

### RENAN sub genome
g=ggplot(d_genome[!is.na(d_genome$fullname),], aes(x=genome, y=percent_DNA, fill=fullname))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette="Paired")
plot(g)

##TE_cumul_length_secale_Lo7.png
### RENAN sub genome
g2=ggplot(d_genome[!is.na(d_genome$fullname),], aes(x=genome, y=cumul_length, fill=fullname))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette="Paired")
plot(g2)




