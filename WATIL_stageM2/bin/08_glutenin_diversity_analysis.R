
library(data.table)
library(tidyverse)

setwd('Y:/ANALYSES/data/storage_proteins_chr1/WATIL_Nezha_stageM2')
dirout='Y:/ANALYSES/results/storage_proteins_chr1/WATIL_Nezha_stageM2/'

seq_length=read_tsv("sequences_length.tsv", col_names = FALSE)
#cons_seq_length = read_tsv("consensus_sequences_length.tsv", col_names = FALSE)
colnames(seq_length)=c("name","length")

cds_coord=read_tsv("18genes_consensus_vs_TaeRenan_coding_part_coord.txt", col_names = F)
colnames(cds_coord)=c("name","start","stop")

vcf_vector = list.files("Y:/ANALYSES/data/STORAGE_PROT_CHR1/WATIL_Nezha_stageM2", pattern='_aling.fasta.vcf')

###debut de boucle sur les fichiers vcf
for (f in 1:length(vcf_vector)) {
  vcf = fread(vcf_vector[f])
  
  Tvcf=t(as.matrix(vcf))
  Tvcf=Tvcf[c(-1,-3,-6,-7,-8,-9),]
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
  Tvcf=gsub(pattern = "/[0-9]:1",replacement = "",Tvcf)
  freq_all=apply(Tvcf[c(-1,-2),-1], MARGIN=2, table)
  POS=names(freq_all)
  
  #recupe nom gene en cours
  nom=sub('\\_aling.fasta.vcf$', '',vcf_vector[f])
  nom=str_replace(nom, 'pbaa_','')
  
  length_amplicon=as.numeric(seq_length[seq_length$name==nom,2])
  

  ###########################################################
  #CALCUL HAPLOTYPES/amplicon entier
  ###########################################################
  
  Tvcf_nuc=as_tibble(Tvcf)
  for (i in 2:dim(Tvcf_nuc)[2]) {
    
    #remplacement allele maj
    Tvcf_nuc[Tvcf_nuc[,i]=="0",i]=Tvcf_nuc[1,i]
    
    #remplacement allele min
    #recuperation de la liste des alleles minoritaires possibles pour la position en cours
    mut_min_vec=unlist(strsplit(as.character(Tvcf_nuc[2,i]),","))
   
    for (k in 1:length(mut_min_vec)) {
      Tvcf_nuc[Tvcf_nuc[,i]==k,i]=mut_min_vec[k]
    }
  }
  #mise en forme: supression lignes 1 et 2, concatenation(unite) des colonnes par ligne
  Tvcf_nuc=Tvcf_nuc[c(-1,-2),]%>%unite("haplo", 2:ncol(Tvcf_nuc), sep="")
  stat_halpo=as.data.frame(table(Tvcf_nuc$haplo))
  stat_halpo$name=paste("haplo",row.names(stat_halpo),sep="_")
  
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
  } else {print(c("no nuc var", nom, f))
    Tvcf_cds=data.frame(Tvcf[,1])}
  
  Tvcf_nuc_cds=as_tibble(Tvcf_cds)
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
    
    ##############################################################
    #EXPORT haplotypes
    write_tsv(data.frame(stat_halpo_cds[,c(3,2,1)]), paste(dirout,nom,"_CDS_",length_cds,"bp_haplotypes.tab", sep=""))
  
    ##############################################################
    #CALCUL Nb SNP/INDEL per 1kb pour l'amplicon CDS
    ##############################################################
    
    nbSNP_per1kb_cds=c(nom, length_cds, ncol(Tvcf_cds)-1, round((ncol(Tvcf_cds)-1)/length_cds*1000, 2))
  
  
    if (f==1) { table_nbSNP_per1kb_cds=nbSNP_per1kb_cds }
    if (f!=1) { table_nbSNP_per1kb_cds=rbind(table_nbSNP_per1kb_cds, nbSNP_per1kb_cds)
    colnames(table_nbSNP_per1kb_cds)=c("name","cds_length", "cds_mut_number","cds_perKb_mut_rate")}
    
  }
  
  
}    ####fin boucle sur les fichiers

table_nbSNP_perKb=left_join(as_tibble(table_nbSNP_per1kb),as_tibble(table_nbSNP_per1kb_cds), by="name")

write_tsv(table_nbSNP_perKb,paste(dirout,"table_nbSNP_perKb.tsv",sep=""))

