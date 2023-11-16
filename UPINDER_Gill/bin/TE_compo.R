library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

#localisation donnees entree
dirin<-"Y:/ANALYSES/bin/soutien_bioinfo/UPINDER_Gill/results/"
setwd(dirin)

classi=read_tsv("Y:/ANALYSES/results/wheatomics_wp1/annotation/TE/TE_classification.txt", col_names=T)

d=read_tsv("Length_TE_AeU_SUPERFamLevel.tsv", col_names=F)
names(d)<-c("chrom", "fam", "length")
ls.str(d)

d=d%>%mutate(chrom=gsub("[Cc]hr", "", chrom), genome=str_extract(chrom, "[:upper:]"))
d=left_join(d, classi, by="fam")


total_dna=4281945474
d_genome=d%>%group_by(fam,fullname)%>%summarise(cumul_length=sum(length), percent_DNA=round(sum(length)/total_dna*100,2))
d_genome=d_genome%>%cbind(genome="U")

##TE_cumul_length_AeU.png
# g2=ggplot(d_genome[!is.na(d_genome$fullname),], aes(x=genome, y=cumul_length, fill=fullname))+
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_brewer(palette="Paired")
# plot(g2)

##################################################################################################
###  FAM level comparison
##################################################################################################

d1=read_tsv("Length_TE_AeU_FamLevel.tsv", col_names=F)
names(d1)<-c("chrom", "fam", "length")
ls.str(d1)
d1=d1%>%mutate(chrom=gsub("[Cc]hr", "", chrom), genome=str_extract(chrom, "[:upper:]"))
dU=d1%>%group_by(fam)%>%summarise(cum_length_U=sum(length), percent_DNA_U=sum(length)/total_dna*100)


d2=read_tsv("Length_TE_refseqv2_FamLevel.tsv", col_names=F)
names(d2)<-c("chrom", "fam", "length")
ls.str(d2)
d2=d2%>%mutate(chrom=gsub("[Cc]hr", "", chrom), genome=str_extract(chrom, "[:upper:]"))

A_dna=4975319984
dA=d2%>%filter(genome=="A")%>%group_by(fam)%>%summarise(cum_length_A=sum(length), percent_DNA_A=sum(length)/A_dna*100)

B_dna=5249122710
dB=d2%>%filter(genome=="B")%>%group_by(fam)%>%summarise(cum_length_B=sum(length), percent_DNA_B=sum(length)/B_dna*100)

D_dna=4001386677
dD=d2%>%filter(genome=="D")%>%group_by(fam)%>%summarise(cum_length_D=sum(length), percent_DNA_D=sum(length)/D_dna*100)

dUA=inner_join(dU, dA, by="fam")%>%mutate(log2FC=log2(percent_DNA_U/percent_DNA_A))%>%
  arrange(log2FC)%>%filter(log2FC>2 | log2FC<(-2))%>%filter(cum_length_U>=100000 | cum_length_A>=100000)

dUB=inner_join(dU, dB, by="fam")%>%mutate(log2FC=log2(percent_DNA_U/percent_DNA_B))%>%
  arrange(log2FC)%>%filter(log2FC>2 | log2FC<(-2))%>%filter(cum_length_U>=100000 | cum_length_B>=100000)

dUD=inner_join(dU, dD, by="fam")%>%mutate(log2FC=log2(percent_DNA_U/percent_DNA_D))%>%
  arrange(log2FC)%>%filter(log2FC>2 | log2FC<(-2))%>%filter(cum_length_U>=100000 | cum_length_D>=100000)

write.table(dUA, file="TEfam_logFC_AeU_versus_CS_A.tsv", sep="\t", col.names=T, row.names=F)
write.table(dUB, file="TEfam_logFC_AeU_versus_CS_B.tsv", sep="\t", col.names=T, row.names=F)
write.table(dUD, file="TEfam_logFC_AeU_versus_CS_D.tsv", sep="\t", col.names=T, row.names=F)
