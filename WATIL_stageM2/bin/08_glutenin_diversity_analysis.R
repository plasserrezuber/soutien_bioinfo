
library(data.table)
library(tidyverse)

seq_lenght <- read.csv("sequences_lenght.csv",sep = "",header = FALSE)
seq_lenght$V1= sub('\\_cons.fasta$', '',seq_lenght$V1)

cds_cord_snp <- fread("18genes_consensus_coding_part_coord.txt")
colnames(cds_cord_snp)=c("gene", "start", "stop")

vcf_vector <- list.files("C:/Users/nwatil/Documents", pattern='_aling.fasta.vcf')

###debut de boucle sur les fichiers vcf
for (f in 1:length(vcf_vector)) {
  vcf <- fread(vcf_vector[f])
  
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
  test=Tvcf%>%mutate(across(starts_with("pos"), ~str_replace(.,"/\\d:\\d","")))
  
  
  
  Tvcf=as.matrix(Tvcf)
  Tvcf=gsub(pattern = "/[0-9]:1",replacement = "",Tvcf)
  freq_all=apply(Tvcf[c(-1,-2),-1], MARGIN=2, table)
  POS=names(freq_all)
  
  #recupe nom gene en cours
  nom=sub('\\_aling.fasta.vcf$', '',vcf_vector[f])
  nom=str_replace(nom, 'pbaa_','') 
  nom=str_replace(nom, 'hmw_','hmw') 
  length=seq_lenght[seq_lenght$V1==nom,2]
  
  ##partie CDS
  
  start=as.numeric(cds_cord_snp[gene==nom,"start"])
  stop=as.numeric(cds_cord_snp[gene==nom,"stop"])
  length_cds=stop-start
  pos_list=data.frame(vcf[vcf$POS>=start & vcf$POS<=stop,"POS"])
  if (dim(pos_list)[1]!=0) {
    pos_list=paste("pos_",pos_list$POS,sep="")
    #Tvcf partie codante
    Tvcf_cds=Tvcf[,c("sample",pos_list)]
  } else   {Tvcf_cds=data.frame(Tvcf[,1])}
  
  ###########################################################
  #CALCUL HAPLOTYPES/amplicon entier
  ###########################################################
  
  Tvcf_nuc=as_tibble(Tvcf)
  for (i in 2:dim(Tvcf_nuc)[2]) {
    
    #remplacement all maj
    Tvcf_nuc[Tvcf_nuc[,i]=="0",i]=Tvcf_nuc[1,i]
    
    #remplacement all min
    mut_min_vec=unlist(strsplit(as.character(Tvcf_nuc[2,i]),","))
   
    for (k in 1:length(mut_min_vec)) {
      Tvcf_nuc[Tvcf_nuc[,i]==k,i]=mut_min_vec[k]
    }
  }
  Tvcf_nuc=Tvcf_nuc[c(-1,-2),]%>%unite("haplo", 2:ncol(Tvcf_nuc), sep="")
  stat_halpo=as.data.frame(table(Tvcf_nuc$haplo))
  stat_halpo$name=paste("haplo",row.names(stat_halpo),sep="_")
  
  ###########################################################
  #CALCUL HAPLOTYPES/région CDS
  ###########################################################
  
  Tvcf_nuc_cds=as_tibble(Tvcf_cds)
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
  #CALCUL Nb SNP per 1kb pour l'amplicon entier et CDS
  ##############################################################
  
  ##Calcul de nombre de mutations sur cds
  nbSNP_per1kb_cds=c(nom, length_cds, ncol(Tvcf_cds)-1, round((ncol(Tvcf_cds)-1)/length_cds*1000, 2))
  
  
  if (f==1) {table_nbSNP_per1kb_cds=nbSNP_per1kb_cds}
  if (f!=1) {table_nbSNP_per1kb_cds=rbind(table_nbSNP_per1kb_cds, nbSNP_per1kb_cds) }
  
  
  ##Calcul de nombre de mutations sur la séquence
  nbSNP_per1kb=c(nom, length, ncol(Tvcf)-1, round((ncol(Tvcf)-1)/length*1000, 2))
  
  
  if (f==1) {table_nbSNP_per1kb=nbSNP_per1kb}
  if (f!=1) {table_nbSNP_per1kb=rbind(table_nbSNP_per1kb, nbSNP_per1kb) }

    
    
  ##############################################################
  #EXPORT VCF
  extention <-paste(nom,".csv",sep="")
  write.table(stat,extention,sep= "\t")
  
}
    ##### fin test if freq_all=list ou matrix
    # nrow(stat)/

    ####fin boucle sur les fichiers
write.table(table_nbSNP_per1kb,"table_nbSNP_per1kb.tsv",sep= "\t")
write.table(table_nbSNP_per1kb_cds,"table_nbSNP_per1kb_cds.tsv",sep= "\t")
