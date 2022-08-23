#!/bin/bash
module load bedtools/2.27.1

# recuperation de la sequence des genes aux bornes des regions d'etude 1AS, 1BS et 1DS chez chinese_spring REFSEQV1
bedtools getfasta -name+ -fi /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-bed /home/palasser/soutien_bioinfo/extraction_3B_REFSEQV1_ClaireGUERIN_210311.bed \
-fo /home/palasser/soutien_bioinfo/extraction_3B_REFSEQV1_ClaireGUERIN_210311.fasta

## extraction 2000bp avant une liste de genes NAC
## changer End of line du fichier texte d'abord (CRLF -> LF)
## 
join -t$'\t' -1 1 -2 4 <(sed 's/01G/02G/' /home/palasser/soutien_bioinfo/488_NAC_CS_RefSeqv1_ID.txt |sed -E 's/\.[0-9]//' |sort) \
<(sort -k4,4 /home/palasser/data/IWGSC_v1.1_20170706.bed) |gawk -v OFS='\t' '{print $2,$3,$3+1,$1}' \
> /home/palasser/soutien_bioinfo/487_NAC_CS_RefSeqv1_start_start.bed

bedtools slop -l 2000 -r 2000 -i /home/palasser/soutien_bioinfo/487_NAC_CS_RefSeqv1_start_start.bed \
-g <(cut -f1,2 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta.fai) \
> /home/palasser/soutien_bioinfo/487_NAC_CS_RefSeqv1_2k_before_2k_afterStart.bed

bedtools getfasta -name -fi /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-bed /home/palasser/soutien_bioinfo/487_NAC_CS_RefSeqv1_2k_before_2k_afterStart.bed \
-fo /home/palasser/soutien_bioinfo/487_NAC_CS_RefSeqv1_2k_before_2k_afterStart.fasta