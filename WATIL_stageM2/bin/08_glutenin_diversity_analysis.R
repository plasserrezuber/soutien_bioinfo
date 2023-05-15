
library(data.table)
library(tidyverse)
library(plyr)
library(broom)

setwd('Y:/ANALYSES/data/storage_proteins_chr1/WATIL_Nezha_stageM2')
dirout='Y:/ANALYSES/results/storage_proteins_chr1/WATIL_Nezha_stageM2/'

seq_length=read_tsv("sequences_length.tsv", col_names = FALSE)
#cons_seq_length = read_tsv("consensus_sequences_length.tsv", col_names = FALSE)
colnames(seq_length)=c("name","length")

cds_coord=read_tsv("18genes_consensus_vs_TaeRenan_coding_part_coord.txt", col_names = F)
colnames(cds_coord)=c("name","start","stop")

indiv=read_tsv("list_indiv_Exige_Glu_PacBio.txt", col_names = F)
list_indiv=rbind("REF","ALT",list_indiv)
colnames(list_indiv)="sample"

vcf_vector = list.files("Y:/ANALYSES/data/storage_proteins_chr1/WATIL_Nezha_stageM2", pattern='_aling.fasta.vcf')

###NB: sans option de consensus, l'outil msa2vcf qui génère les vcf à partir d'un alignement multiple considere comme
# chrUn la longueur totale de l'alignement multiple
# les coordonnees sont celles sur le chrUn
# l'allele REF est l'allele majoritaire (="0"), les alleles ALT (1, 2, 3, etc) sont les alleles minoritaires

###debut de boucle sur les fichiers vcf
for (f in 1:length(vcf_vector)) {
  vcf = fread(vcf_vector[f])
#   vcf=vcf%>%
#     select("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",contains("cluster-0_ReadCount"))%>%
#     rename_at(vars(starts_with("sample")), funs(str_sub(.,8,13)))
#   
#   #recupe nom gene en cours et longueur amplicon en cours
#   nom=sub('\\_aling.fasta.vcf$', '',vcf_vector[f])
#   nom=str_replace(nom, 'pbaa_','')
#   length_amplicon=as.numeric(seq_length[seq_length$name==nom,2])
#   
#   vcf$POS=paste(nom,vcf$POS, sep="_")
#   
#   #liste de colonne pour concatener les vcf
#   listcol=as.vector(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",as_vector(indiv)))
#   tibcol=as_tibble(matrix(nrow=1,ncol=length(listcol), dimnames=list(1,listcol)))
#   
#   if (f==1) {table_vcf=rbind.fill(tibcol, vcf)[-1,]}
#   if (f!=1) {table_vcf=rbind.fill(table_vcf, vcf) }
# 
#   fwrite(table_vcf,"PacBio_94indiv_18_glutenin_Exige.vcf", sep="\t")
# }
  
  vcf=vcf%>%select("POS","REF","ALT",contains("cluster-0_ReadCount"))
  
  Tvcf=t(as.matrix(vcf))
  Tvcf=cbind(rownames(Tvcf),Tvcf)
  Tvcf=data.frame(Tvcf)
  
  colnames(Tvcf)=gsub(pattern = " ",replacement = "",Tvcf[1,])
  colnames(Tvcf)=paste("pos_",colnames(Tvcf),sep="")
  colnames(Tvcf)[1]="sample"
  Tvcf=Tvcf[-1,]
  rownames(Tvcf)=NULL
  Tvcf=tibble(Tvcf)
  Tvcf=Tvcf%>%mutate(sample=str_extract(sample,"bc20.."))
  Tvcf[1,1]="REF"
  Tvcf[2,1]="ALT"
  #Tvcf=Tvcf%>%mutate(across(starts_with("pos"), ~str_replace(.,"/\\d:\\d","")))
  
  Tvcf=as.matrix(Tvcf)
  Tvcf=gsub(pattern = "/[0-9]*:1",replacement = "",Tvcf)
  Tvcf=gsub(pattern = "[.]/[.]",replacement = NA, Tvcf)
  

  ##############################################################
  ### traitement pour export table complete (code 0/1/2) pour les 18 genes codant les glutenines
  
  ## pour etre sur que la liste des indiv soit les memes pour chaque gene et faire un cbind ensuite.
  Tvcf=as_tibble(Tvcf)
  Tvcf_GENO=left_join(list_indiv,Tvcf, by="sample")
  
  colnames(Tvcf_GENO)=gsub(pattern="pos", replacement=nom,colnames(Tvcf_GENO))
  
  if (f==1) {table_GENO=Tvcf_GENO}
  if (f!=1) {table_GENO=cbind(table_GENO, Tvcf_GENO[,-1])}
  
  
  ##############################################################
  ## Calcul des frequences alleliques et pourcentage NA
  
  freq_allelic=Tvcf_GENO%>%filter(sample!="REF" & sample!="ALT")%>%select(-sample)%>%
    pivot_longer(everything(), names_to = "pos", values_to="allele")%>%
    group_by(pos,allele)%>%
    summarise(n=n())%>%
    mutate(freq_perc=n/nrow(Tvcf_GENO)*100)%>%
    mutate(pos=str_replace(pos, "pos", nom))
  
  MAF_NA=left_join(freq_allelic%>%group_by(pos)%>%summarise(MAF=min(freq_perc)),
                   freq_allelic%>%filter(is.na(allele))%>%select(pos,freq_perc)%>%rename(NAperc=freq_perc), by="pos")
  
  if (f==1) {table_MAF_NA=MAF_NA}
  if (f!=1) {table_MAF_NA=rbind(table_MAF_NA, MAF_NA)}

  ###########################################################
  #CALCUL HAPLOTYPES/amplicon entier
  ###########################################################
  
  ## transformation matrice geno codage chiffres en codage nucleotides
  
  Tvcf_nuc=Tvcf
  for (i in 2:dim(Tvcf_nuc)[2]) {
    
    #remplacement allele maj=REF
    Tvcf_nuc[Tvcf_nuc[,i]=="0" & !is.na(Tvcf_nuc[,i]),i]=Tvcf_nuc[1,i]
    
    #remplacement alleles min=ALT par 
    #recuperation de la liste des alleles minoritaires possibles pour la position en cours
    mut_min_vec=unlist(strsplit(as.character(Tvcf_nuc[2,i]),","))
   
    for (k in 1:length(mut_min_vec)) {
      Tvcf_nuc[Tvcf_nuc[,i]==k & !is.na(Tvcf_nuc[,i]),i]=mut_min_vec[k]
    }
  }
  
  ##############################################################
  ### traitement pour export table complete (code ATGC) pour les 18 genes codant les glutenines
  
  ## pour etre sur que la liste des indiv soit les memes pour chaque gene et faire un cbind ensuite.
  Tvcf_GENO_nuc=left_join(list_indiv,Tvcf_nuc, by="sample")
  
  colnames(Tvcf_GENO_nuc)=gsub(pattern="pos", replacement=nom,colnames(Tvcf_GENO_nuc))
  
  if (f==1) {table_GENO_nuc=Tvcf_GENO_nuc}
  if (f!=1) {table_GENO_nuc=cbind(table_GENO_nuc, Tvcf_GENO_nuc[,-1])}

  
  
  #mise en forme: supression lignes 1 et 2, concatenation(unite) des colonnes par ligne
  Tvcf_nuc=Tvcf_nuc[c(-1,-2),]%>%unite("haplo", 2:ncol(Tvcf_nuc), sep="")
  stat_halpo=as.data.frame(table(Tvcf_nuc$haplo))
  stat_halpo$name=paste("haplo",row.names(stat_halpo),sep="_")
  
  print(nom)
  print(sum(stat_halpo$Freq))
  
  ##############################################################
  #EXPORT haplotypes
  write_tsv(data.frame(stat_halpo[,c(3,2,1)]), paste(dirout,nom,"_amplicon_",length_amplicon,"bp_haplotypes.tab", sep=""))
  
  
  ##############################################################
  #CALCUL Nb SNP/INDEL per 1kb pour l'amplicon entier
  ##############################################################
  
  nbSNP_per1kb=c(nom, length_amplicon, ncol(Tvcf)-1, round((ncol(Tvcf)-1)/length_amplicon*1000, 2))
  
  
  if (f==1) {table_nbSNP_per1kb=nbSNP_per1kb}
  if (f!=1) {table_nbSNP_per1kb=rbind(table_nbSNP_per1kb, nbSNP_per1kb)
  colnames(table_nbSNP_per1kb)=c("name","amplicon_length", "mut_number","perKb_mut_rate")}
  
  ###########################################################
  #CALCUL HAPLOTYPES/region CDS
  ###########################################################
  
  
  ####PB revoir script gmap 18 cds ref contre gmapdb consensus 18 amplicons assembles pacbio
  start=as.numeric(cds_coord[cds_coord$name==nom,"start"])
  stop=as.numeric(cds_coord[cds_coord$name==nom,"stop"])
  length_cds=stop-start
  pos_list=data.frame(vcf[vcf$POS>=start & vcf$POS<=stop,"POS"])
  
  if (dim(pos_list)[1]!=0) {
    pos_list=paste("pos_",pos_list$POS,sep="")
    #Tvcf partie codante
    Tvcf_cds=Tvcf[,c("sample",pos_list)]
  } else {print(c("no nucleotidic variability in CDS", nom, f))
    Tvcf_cds=data.frame(Tvcf[,1])
    nbSNP_per1kb_cds=c(nom, length_cds, ncol(Tvcf_cds)-1, (ncol(Tvcf_cds)-1)/length_cds*1000)
    }
  
  Tvcf_nuc_cds=as_tibble(Tvcf_cds)
  
  
  # if il y a au moins un polymorphisme dans la cds
  if (dim(Tvcf_nuc_cds)[2]>=2) {
    for (i in 2:dim(Tvcf_nuc_cds)[2]) {
    
      #remplacement all maj
      Tvcf_nuc_cds[Tvcf_nuc_cds[,i]=="0",i]=Tvcf_nuc_cds[1,i]
    
      #remplacement all min
      mut_min_cds=unlist(strsplit(as.character(Tvcf_nuc_cds[2,i]),","))
    
    
      for (k in 1:length(mut_min_cds)) {
        Tvcf_nuc_cds[Tvcf_nuc_cds[,i]==k,i]=mut_min_cds[k]
      }
    }
    Tvcf_nuc_cds=Tvcf_nuc_cds[c(-1,-2),]%>%unite("haplo", 2:ncol(Tvcf_nuc_cds), sep="")
    stat_halpo_cds=as.data.frame(table(Tvcf_nuc_cds$haplo))
    stat_halpo_cds$name=paste("haplo",row.names(stat_halpo_cds),sep="_")
    
    print(nom)
    print(sum(stat_halpo$Freq))
    
    ##############################################################
    #EXPORT haplotypes
    write_tsv(data.frame(stat_halpo_cds[,c(3,2,1)]), paste(dirout,nom,"_CDS_",length_cds,"bp_haplotypes.tab", sep=""))
  
    ##############################################################
    #CALCUL Nb SNP/INDEL per 1kb pour l'amplicon CDS
    ##############################################################
    
    nbSNP_per1kb_cds=c(nom, length_cds, ncol(Tvcf_cds)-1, round((ncol(Tvcf_cds)-1)/length_cds*1000, 2))
  
  } # fin de if il y a au moins un polymorphisme
  
  
  if (f==1) { table_nbSNP_per1kb_cds=nbSNP_per1kb_cds }
  if (f!=1) { table_nbSNP_per1kb_cds=rbind(table_nbSNP_per1kb_cds, nbSNP_per1kb_cds)}
  
  
}    ####fin boucle sur les fichiers

colnames(table_nbSNP_per1kb_cds)=c("name","cds_length", "cds_mut_number","cds_perKb_mut_rate")
table_nbSNP_perKb=left_join(as_tibble(table_nbSNP_per1kb),as_tibble(table_nbSNP_per1kb_cds), by="name")

write_tsv(table_nbSNP_perKb,paste(dirout,"table_nbSNP_perKb.tsv",sep=""))

write_tsv(table_GENO,"Y:/ANALYSES/data/storage_proteins_chr1/Exige_FSOV/GENO_PacBio_94indiv_18_glutenin_Exige.tsv")
write_tsv(table_GENO_nuc,"Y:/ANALYSES/data/storage_proteins_chr1/Exige_FSOV/GENO_nuc_PacBio_94indiv_18_glutenin_Exige.tsv")
