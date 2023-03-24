#!/bin/bash
module load bedtools/2.27.1

OUTPUT='/home/palasser/soutien_bioinfo/BAZILE'
mkdir $OUTPUT

## DEMANDE: Recup des sequences des genes RECQ4 de REFSEQV1 suivant chez Renan: TraesCS2A02G304900 TraesCS2B02G321700 TraesCS2D02G303500

###########################################################################################################################################
#### REFSEQV2
###########################################################################################################################################
## recup des geneID REFSEQV2
egrep 'TraesCS2A02G304900'$'\t|TraesCS2B02G321700'$'\t|TraesCS2D02G303500'$'\t' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_IDmaping.csv \
> $OUTPUT/RECQ4_chr2_REFSEQV2_geneID.txt

## VERIF annotation fonctionnelle
#egrep 'TraesCS2A01G304900.1|TraesCS2B01G321700.1|TraesCS2D01G303500.1' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.0/FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB


## VERIF si gene HC et qualite du transfert magatt de REFSEQV1 vers REFSEQV2
liste=$(cut -f2 $OUTPUT/RECQ4_chr2_REFSEQV2_geneID.txt |tr '\n' '|' |sed 's/|$//')
egrep $liste /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC.gff3 |fgrep -w gene \
> $OUTPUT/RECQ4_chr2_REFSEQV2.gff3
## gff3 to bed
gawk -v OFS='\t' '{print $1,$4-1,$5,$9}' $OUTPUT/RECQ4_chr2_REFSEQV2.gff3 |cut -d';' -f1 > $OUTPUT/RECQ4_chr2_REFSEQV2_gene.bed
## bed REFSEQV1
egrep 'TraesCS2A02G304900$|TraesCS2B02G321700$|TraesCS2D02G303500$' ~/data/IWGSC_v1.1_20170706.bed > $OUTPUT/RECQ4_chr2_REFSEQV1_geneID.bed

#extraction fasta gene REFSEQV2
bedtools getfasta -name+ -fi /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/CS_pesudo_v2.1.fa \
-bed $OUTPUT/RECQ4_chr2_REFSEQV2_gene.bed -fo $OUTPUT/RECQ4_chr2_IWGSC_refseqv2.1_annotation_200916_HC_gene.fasta


###########################################################################################################################################
## REFSEQV2 fasta
rm $OUTPUT/RECQ4_chr2_IWGSC_refseqv2.1_annotation_200916_HC_cds.fasta $OUTPUT/RECQ4_chr2_IWGSC_refseqv2.1_annotation_200916_HC_pep.fasta
for g in $(cut -f2 $OUTPUT/RECQ4_chr2_REFSEQV2_geneID.txt);
do
    gawk -v gene=$g 'BEGIN { RS=">" } { if ($0~gene) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC_cds.fasta \
    >> $OUTPUT/RECQ4_chr2_IWGSC_refseqv2.1_annotation_200916_HC_cds.fasta

    gawk -v gene=$g 'BEGIN { RS=">" } { if ($0~gene) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC_pep.fasta \
    >> $OUTPUT/RECQ4_chr2_IWGSC_refseqv2.1_annotation_200916_HC_pep.fasta
done

###########################################################################################################################################
#### RENAN CDS PEP
###########################################################################################################################################
## recup des geneID RENAN FINAUX
egrep $liste /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/TaeRenan_refseq_v2.0_IWGSCv2.1_correspondance.txt \
> $OUTPUT/RECQ4_chr2_RENAN_geneID.txt

## 
listeR=$(cut -f1 $OUTPUT/RECQ4_chr2_RENAN_geneID.txt |tr '\n' '|' |sed 's/|$//')
egrep $listeR /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/TaeRenan_refseq_v2.0_genesHC.gtf \
> $OUTPUT/RECQ4_chr2_RENAN.gtf

rm $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_genesHC_cds.fasta $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_genesHC_pep.fasta

for g in $(cut -f1 $OUTPUT/RECQ4_chr2_RENAN_geneID.txt);
do
    gawk -v gene=$g 'BEGIN { RS=">" } { if ($0~gene) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/TaeRenan_refseq_v2.0_genesHC_cds.fa \
    >> $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_genesHC_cds.fasta

    gawk -v gene=$g 'BEGIN { RS=">" } { if ($0~gene) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/TaeRenan_refseq_v2.0_genesHC_pep.fa \
    >> $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_genesHC_pep.fasta
done



###########################################################################################################################################
#### RENAN GENE
###########################################################################################################################################
egrep $liste /home/herimbert/work/wheatomics/wp1-renan/annot/pseudo_v2/magatt-master/TaeRenan_magatt_2104/TaeRenan_magatt_2104_HC.gff3 |fgrep -w gene \
> $OUTPUT/RECQ4_chr2_RENAN_gene.gff3
## gff3 to bed
gawk -v OFS='\t' '{print $1,$4-1,$5,$9}' $OUTPUT/RECQ4_chr2_RENAN_gene.gff3 |cut -d';' -f1 |sed 's/ID=/previous_ID=/' \
|paste --delimiters=' ' - <(cut -f1 $OUTPUT/RECQ4_chr2_RENAN_geneID.txt) $OUTPUT/RECQ4_chr2_REFSEQV2_geneID.txt \
|tr ' ' '\t' |gawk -v OFS='\t' '{print $1,$2,$3,"ID="$5";"$4,"REFSEQV2_ID="$7";REFSEQV1_ID="$6}' \
> $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_gene.bed


#extraction fasta gene RENAN
bedtools getfasta -name+ -fi /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/assembly/pseudo_v2/Renan_v13_v2.pseudo.v2.fa \
-bed $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_gene.bed -fo $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_genesHC_gene.fasta

sed -i 's/ID=/previous_ID=/' $OUTPUT/RECQ4_chr2_TaeRenan_refseq_v2.0_genesHC_gene.fasta
