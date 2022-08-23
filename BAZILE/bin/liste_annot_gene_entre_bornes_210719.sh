#!/bin/bash
module load bedtools/2.27.1

OUTPUT='/home/palasser/soutien_bioinfo/BAZILE'
mkdir $OUTPUT

## DEMANDE: Extraction liste de genes chez Renan_pseudov2 entre les bornes suivantes de REFSEQV1 (fichier $OUTPUT/ROI_CS_210719.txt):
# chr2A : 514748482..530569798
# chr2B : 454341715..468445368
# chr2D : 375144056..395621869

# fichier: $OUTPUT/ROI_Renan.txt
#chr2A   514748482   530569798
#chr2B   454341715   468445368
#chr2D   375144056   395621869


##### alignement bwa des isbp refseqv1 sur Renan
# sbatch /home/palasser/bin/bwa_isbp_REFSEQV1_vs_RENAN_PSEUDOV2.sh

join -t$'\t' -1 4 -2 4 <(sort -k4,4 /home/palasser/data/REFSEQV1/isbp_refseqv1.bed) \
<(sort -k4,4 /home/palasser/results/bwa/isbp_REFSEQV1_vs_RENAN_PSEUDOV2/ISBPS_REFSEQV1_vs_RENAN_PSEUDOV2_filtered.bed |cut -f1-4) \
|sort -V -k2,2 -k3,3 > /home/palasser/data/RENAN_v2_pseudo/coord_isbp_refseqv1_TaeRenan_v2.0.txt


#passage coord REFSEQV1 a RENAN
rm $OUTPUT/ROI_Renan_210719.txt
while read line;
do
    #coord refseqv1
    chrom=$(echo $line |cut -d' ' -f1 |sed 's/chr//')
    echo $chrom
    start=$(echo $line |cut -d' ' -f2)
    echo $start
    end=$(echo $line |cut -d' ' -f3)
    echo $end
    #coord refseqv2
    start2=$(gawk -v OFS='\t' -v chr=$chrom -v strt=$start -v nd=$end '{ if ($2~chr && $5~chr && $3>=strt && $4<=nd) print $0 }' /home/palasser/data/RENAN_v2_pseudo/coord_isbp_refseqv1_TaeRenan_v2.0.txt |head -n1 |cut -f6)
    echo $start2
    end2=$(gawk -v OFS='\t' -v chr=$chrom -v strt=$start -v nd=$end '{ if ($2~chr && $5~chr && $3>=strt && $4<=nd) print $0 }' /home/palasser/data/RENAN_v2_pseudo/coord_isbp_refseqv1_TaeRenan_v2.0.txt |tail -n1 |cut -f7)
    echo $end2
    #coord TaeRenan_refseq_v2.0
    echo -e "chr${chrom}\t${start2}\t${end2}" >> $OUTPUT/ROI_Renan_210719.txt 
done < $OUTPUT/ROI_CS_210719.txt 

###########################################################################################################################################
#### Extraction info
###########################################################################################################################################
## merge gff HC et LC de Renan:
ml GenomeTools
gt gff3 -sort -tidy -retainids /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/v2021-09-08/TaeRenan_refseqv2.0_genesHC.gff3 \
/storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/v2021-09-08/TaeRenan_refseqv2.0_genesLC.gff3 \
1>  /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_genes.gff3 2> /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_genes_gt.log


rm $OUTPUT/TaeRenan_refseq_v2.0_genes_ROI_RECQ4.gff3
while read line;
do
    chrom=$(echo $line |cut -d' ' -f1 |sed 's/chr//')
    start=$(echo $line |cut -d' ' -f2)
    end=$(echo $line |cut -d' ' -f3)
    echo -e "##sequence-region\t$line" >> $OUTPUT/TaeRenan_refseq_v2.0_genes_ROI_RECQ4.gff3
    gawk -v OFS='\t' -v chr=$chrom -v strt=$start -v nd=$end '{ if ($1~chr && $4>=strt && $5<=nd) print $0 }' /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseq_v2.0_genes.gff3 \
    >> $OUTPUT/TaeRenan_refseq_v2.0_genes_ROI_RECQ4.gff3
done < $OUTPUT/ROI_Renan_210719.txt 

fgrep -w 'gene' $OUTPUT/TaeRenan_refseq_v2.0_genes_ROI_RECQ4.gff3 |gawk -v FS='ID=' '{print $2}' |tr ';' '\t' |gawk -v OFS='\t' 'NF>0 {for (i=1;i<=NF;i++) { if ($i~"previous_id") print $1,$i } }' \
|sed 's/previous_id=//' |cut -d',' -f1 > $OUTPUT/liste_TaeRenan_refseq_v2.0_REFSEQV2_geneID_ROI_RECQ4.txt


join -t$'\t' -1 2 -2 2 <(sort -k2,2 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_IDmaping.csv) \
<(sort -k2,2 $OUTPUT/liste_TaeRenan_refseq_v2.0_REFSEQV2_geneID_ROI_RECQ4.txt) \
|sed 's/02G/01G/' |gawk -v OFS='\t' '{print $1,$2".1",$3}' > $OUTPUT/liste_REFSEQV2_REFSEQV1_TaeRenan_refseq_v2.0_geneID_ROI_RECQ4.txt

cat <(join -t$'\t' -1 2 -2 1 <(sort -k2,2 $OUTPUT/liste_REFSEQV2_REFSEQV1_TaeRenan_refseq_v2.0_geneID_ROI_RECQ4.txt) \
<(sort -k1,1 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.0/FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB)) \
<(join -t$'\t' -1 2 -2 1 <(sort -k2,2 $OUTPUT/liste_REFSEQV2_REFSEQV1_TaeRenan_refseq_v2.0_geneID_ROI_RECQ4.txt) \
<(sort -k1,1 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.0/FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__LCgenes_v1.0.TAB)) \
> $OUTPUT/functional_annotation_REFSEQV1_TaeRenan_refseq_v2.0_ROI_RECQ4.txt