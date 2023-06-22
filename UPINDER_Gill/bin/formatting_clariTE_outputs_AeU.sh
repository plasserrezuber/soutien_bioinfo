#!/bin/bash

module load bedtools/2.27.1 GenomeTools

DATA='/home/palasser/projects/soutien_bioinfo/UPINDER_Gill/data'

############################################################################################################################################################################
###### MAKE MERGED GFF FOR WHOLE GENOME
OUTPUT='/home/palasser/projects/soutien_bioinfo/UPINDER_Gill/results/annotTE_AeU'
gt gff3 -sort -tidy -retainids $OUTPUT/AeU_RefSeq_*_clariTE_friendly.gff3 1>  $OUTPUT/AeU_RefSeq_clariTE.gff3 2> $OUTPUT/AeU_RefSeq_clariTE_gt.log

gt gff3validator $OUTPUT/AeU_RefSeq_clariTE_friendly.gff3

###### gff3 to bed
egrep '\srepeat_region' $OUTPUT/AeU_RefSeq_clariTE.gff3 |gawk -v OFS='\t' '{print $1,$4-1,$5,$9}' > $OUTPUT/AeU_RefSeq_clariTE.bed

###### masking before Triannot
# ml bedtools
# bedtools maskfasta -fi $DATA/FINAL_AeU_asm.fasta -bed $OUTPUT/AeU_RefSeq_clariTE.bed -fo $OUTPUT/AeU_RefSeq_TEmasked.fasta

# ml exonerate pour verif taux N dans genome
## genome entier
# fastacomposition $OUTPUT/AeU_RefSeq_TEmasked.fasta
## chrom par chrom
# fastacomposition -s --fasta $OUTPUT/AeU_RefSeq_TEmasked.fasta


gawk '{ if ($0!~"ChrUn") sum+=$2; print sum}' /home/palasser/projects/soutien_bioinfo/UPINDER_Gill/data/FINAL_AeU_asm.fasta.fai |tail -n1
### 4281945474 GENOME LENGTH

###### TOTAL GENOME LENGTH repeat_region:
## bedtools merge: Merges overlapping BED/GFF/VCF entries into a single interval, the input file (-i) file must be sorted by chrom, then start.
bedtools merge -i <(grep '\srepeat_region\s' $OUTPUT/AeU_RefSeq_clariTE.gff3 |grep -v '^ChrUn' |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 3 669 965 436

###### TOTAL GENOME LENGTH match_part:
bedtools merge -i <(grep '\smatch_part\s' $OUTPUT/AeU_RefSeq_clariTE.gff3 |grep -v '^ChrUn' |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 3555472160 (diff avec repeat_region s'explique par le fait que les match_part ne se suivent pas a la base pres (intervalle entre deux par ex), alors que repeat_region englobe tout)

###### TOTAL GENOME LENGTH ALL FEATURES (match_part overlapping repeat_region):
bedtools merge -i <(grep '^Chr' $OUTPUT/AeU_RefSeq_clariTE.gff3 |grep -v '^ChrUn' |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 3669965541 (diff avec repeat_region=105)

/home/palasser/projects/wheatomics_wp1/annotation/TE/bin/length_Family_TE.py $OUTPUT/AeU_RefSeq_clariTE.gff3 |sort -k1,1 -k4,4n > $OUTPUT/AeU_RefSeq_clariTE_For_Length_TE_calculation.gff3

for chr in "1U" "2U" "3U" "4U" "5U" "6U" "7U" "Un";
do
    chrom='Chr'$chr
    # grep ^$chrom $OUTPUT/AeU_RefSeq_clariTE_For_Length_TE_calculation.gff3 |gawk -v OFS='\t' '{print $NF,$5-$4}' |sed -E 's/\.[0-9]*//' \
    # |gawk -v grp=$chrom -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print grp,i,a[i]}'  |sort -k1,1 >> Length_TE_AeU_FamLevel.tsv

    grep ^$chrom $OUTPUT/AeU_RefSeq_clariTE_For_Length_TE_calculation.gff3 |gawk -v OFS='\t' '{print $NF,$5-$4}' |sed -E 's/_.*\t/\t/' \
    |gawk -v grp=$chrom -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print grp,i,a[i]}'  |sort -k1,1 >> Length_TE_AeU_SUPERFamLevel.tsv

    grep ^$chrom $OUTPUT/AeU_RefSeq_clariTE.gff3 |grep '\srepeat_region\s' |cut -f9 |sed -E 's/ID=.*Family:(.*)[_|w].* Matching.*/\1/' \
    |sort |uniq -c |gawk -v grp=$chrom -v OFS='\t' '{print grp,$2,$1}' >> Nb_TE_AeU_SUPERFamLevel.tsv
done

##CALCUL A RETENIR
gawk '{if ($0!~"ChrUn") sum+=$3; print sum}' Length_TE_AeU_SUPERFamLevel.tsv |tail -n1
### 3384737944
gawk '{ if ($0!~"ChrUn" && $0!~"unknown") sum+=$3; print sum}' Length_TE_AeU_SUPERFamLevel.tsv |tail -n1
### 3383953312