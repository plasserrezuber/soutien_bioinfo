#!/bin/bash

# Dans R extrait de 08_glutenin_diversity_analysis.R:
# library(data.table)
# library(tidyverse)
# library(broom)

# setwd('Y:/ANALYSES/data/storage_proteins_chr1/WATIL_Nezha_stageM2')
# dirout='Y:/ANALYSES/results/storage_proteins_chr1/WATIL_Nezha_stageM2/'

# seq_length=read_tsv("sequences_length.tsv", col_names = FALSE)
# #cons_seq_length = read_tsv("consensus_sequences_length.tsv", col_names = FALSE)
# colnames(seq_length)=c("name","length")

# cds_coord=read_tsv("18genes_consensus_vs_TaeRenan_coding_part_coord.txt", col_names = F)
# colnames(cds_coord)=c("name","start","stop")

# list_indiv=read_tsv("list_indiv_Exige_Glu_PacBio.txt", col_names = F)
# list_indiv=rbind("REF","ALT",list_indiv)
# colnames(list_indiv)="sample"

# vcf_vector = list.files("Y:/ANALYSES/data/storage_proteins_chr1/WATIL_Nezha_stageM2", pattern='_aling.fasta.vcf')

# for (f in 1:length(vcf_vector)) {
#   vcf = fread(vcf_vector[f])
#   vcf=vcf%>%
#     select("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",contains("cluster-0_ReadCount"))%>%
#     rename_at(vars(starts_with("sample")), funs(str_sub(.,8,13)))
  
#   #recupe nom gene en cours et longueur amplicon en cours
#   nom=sub('\\_aling.fasta.vcf$', '',vcf_vector[f])
#   nom=str_replace(nom, 'pbaa_','')
#   length_amplicon=as.numeric(seq_length[seq_length$name==nom,2])
  
#   vcf$POS=paste(nom,vcf$POS, sep="_")
  
#   #liste de colonne pour concatener les vcf
#   listcol=as.vector(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",as_vector(indiv)))
#   tibcol=as_tibble(matrix(nrow=1,ncol=length(listcol), dimnames=list(1,listcol)))
  
#   if (f==1) {table_vcf=rbind.fill(tibcol, vcf)[-1,]}
#   if (f!=1) {table_vcf=rbind.fill(table_vcf, vcf) }

#   fwrite(table_vcf,"PacBio_94indiv_18_glutenin_Exige.vcf", sep="\t")
# }


#############################################################################################
### Strategie abandonnee pour finalement concatener les vcf sous R
ml gcc/8.1.0 bcftools/1.9

## recup de fichiers vcf après retrait des colonnes contenant "cluster1_" sous R
## recup des header ##
for file in $(find pbaa_*_aling.fasta.vcf);
do
    f=$(cut -d'.' -f1 <(echo $file))
    grep '##' $file > cluster0_${f}.vcf
    cat cluster0_${f}.vcf cluster0_${f}.fasta.vcf > tmp && mv tmp cluster0_${f}.vcf
done

#pour harmoniser les noms des samples:
for f in $(find cluster0_*_aling.vcf);
do
     sed 's/ //g' $f |sed -E 's/\tsample-(bc20[0-9]{2})[^\t]*/\t\1\t/g' > tmp && mv tmp $f
done

# bcftools viwe pour spliter par sample (par gene)
# bcftools concat par sample tous les genes
# bcftools merge tous les genes
#############################################################################################