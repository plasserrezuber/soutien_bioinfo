#!/bin/bash
module load bedtools/2.27.1

OUTPUT='/home/palasser/soutien_bioinfo/BOUDET'
mkdir $OUTPUT

## DEMANDE: annotation detaillee des 3 genes suivants:
TraesCS7A02G569100
TraesCS7B02G489500
TraesCS7D02G543500

####################################################################################
### recup GFF REFSEQV1
####################################################################################
egrep 'TraesCS7A02G569100|TraesCS7B02G489500|TraesCS7D02G543500' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/IWGSC_v1.1_HC_20170706.gff3) \
> $OUTPUT/IWGSC_v1.1_HC_3genes_BOUDET.gff3

########################################################################################################################################################################
########################################################################################################################################################################
### correspondance REFSEQV2
TraesCS7A03G1388600
TraesCS7B03G1341000
TraesCS7D03G1305700

####################################################################################
### recup GFF REFSEQV2
####################################################################################
egrep 'TraesCS7A03G1388600|TraesCS7B03G1341000|TraesCS7D03G1305700' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC.gff3) \
> $OUTPUT/IWGSC_v2.1_HC_3genes_BOUDET.gff3

####################################################################################
# ml gcc/8.1.0 cufflinks/2.2.1
# for g in "TraesCS7A03G1388600" "TraesCS7B03G1341000" "TraesCS7D03G1305700";
# do
#     gffread <(grep $g $OUTPUT/IWGSC_v2.1_HC_3genes_BOUDET.gff3) \
#     -g /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/CS_pesudo_v2.1.fa \
#     -w $OUTPUT/IWGSC_v2.1_HC_${g}_exon.fasta
# done

# for g in "TraesCS7A03G1388600" "TraesCS7B03G1341000" "TraesCS7D03G1305700";
# do
#     gawk -v seq=$g 'BEGIN { RS=">" } { if ($0 ~ seq) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC_mrna.fasta \
#     > $OUTPUT/IWGSC_v2.1_HC_${g}_mrna.fasta
# done

####################################################################################
#### gff to bed
####################################################################################
gawk -F'\t' '{print $1"\t"$4-1"\t"$5"\t"$9}' <(grep 'gene' $OUTPUT/IWGSC_v2.1_HC_3genes_BOUDET.gff3) \
|cut -d ';' -f1 | sed "s/ID\=//" |sort -V -k1,1 -k2,2 > $OUTPUT/IWGSC_v2.1_HC_3genes_BOUDET.bed

####################################################################################
## selection reads overlapping nac   !!!!!!!!!!!!!!!!!
####################################################################################
## {.} permet de garder la variable introduite par ::: sans l'extension ( avec chemin+nom de fichier)
sbatch -p fast -c 8 --wrap="ml samtools parallel; parallel -j 8 'samtools view -bL /home/palasser/soutien_bioinfo/BOUDET/IWGSC_v2.1_HC_3genes_BOUDET.bed {} > {.}_NAC.bam' ::: /home/palasser/results/hisat2/IWGSC_refseqv2.1_mRNA/*_q60.bam"
sbatch -p fast -c 8 --wrap="ml samtools; samtools merge -@8 /home/palasser/soutien_bioinfo/BOUDET/NAC_merge.bam /home/palasser/results/hisat2/IWGSC_refseqv2.1_mRNA/*q60_NAC.bam"
samtools index -c NAC_merge.bam


#### charger les fichiers suivant dans IGV:
# $OUTPUT/IWGSC_v2.1_HC_3genes_BOUDET.bed
# HMW_merge.bam
# HMW_merge.bam.csi
# /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa
# /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0.fa.fai