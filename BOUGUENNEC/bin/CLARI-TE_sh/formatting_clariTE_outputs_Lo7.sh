#!/bin/bash

module load bedtools/2.27.1 GenomeTools

DATA='/storage/groups/gdec/shared/Secale_cereale/Lo7_RefSeq'  

############################################################################################################################################################################
###### MAKE MERGED GFF FOR WHOLE GENOME
OUTPUT='/home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7'
gt gff3 -sort -tidy -retainids $OUTPUT/Lo7_RefSeq_*_clariTE_friendly.gff3 1>  $OUTPUT/Lo7_RefSeq_clariTE.gff3 2> $OUTPUT/Lo7_RefSeq_clariTE_gt.log

###### gff3 to bed
egrep '\srepeat_region' $OUTPUT/Lo7_RefSeq_clariTE.gff3 |gawk -v OFS='\t' '{print $1,$4-1,$5,$9}' > $OUTPUT/Lo7_RefSeq_clariTE.bed

###### masking before Triannot
# ml bedtools
# bedtools maskfasta -fi $DATA/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta -bed $OUTPUT/Lo7_RefSeq_clariTE.bed -fo $OUTPUT/Lo7_RefSeq_TEmasked.fasta

###### generate TEs fasta file  !!!!VERIFIER SI JUSTE!!!!
# ml bedtools
#bedtools getfasta -fi $DATA/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta -bed $OUTPUT/Lo7_RefSeq_clariTE.bed -fo $OUTPUT/Lo7_RefSeq_clariTE.fasta

# ml exonerate pour verif taux N dans genome
## genome entier
# fastacomposition $OUTPUT/Lo7_RefSeq_TEmasked.fasta
## chrom par chrom
# fastacomposition -s --fasta $OUTPUT/Lo7_RefSeq_TEmasked.fasta


gawk '{ if ($0!~"chrUn") sum+=$2; print sum}' /storage/groups/gdec/shared/Secale_cereale/Lo7_RefSeq/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta.fai
### 6206789216 GENOME LENGTH

###### TOTAL GENOME LENGTH repeat_region:
## bedtools merge: Merges overlapping BED/GFF/VCF entries into a single interval, the input file (-i) file must be sorted by chrom, then start.
bedtools merge -i <(grep '\srepeat_region\s' $OUTPUT/Lo7_RefSeq_clariTE.gff3 |grep -v '^chrUn' |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 5561966586 

###### TOTAL GENOME LENGTH match_part:
bedtools merge -i <(grep '\smatch_part\s' $OUTPUT/Lo7_RefSeq_clariTE.gff3 |grep -v '^chrUn' |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 5376806824 (diff avec repeat_region s'explique par le fait que les match_part ne se suivent pas a la base pres (intervalle entre deux par ex), alors que repeat_region englobe tout)
##CALCUL A RETENIR
gawk '{if ($0!~"chrUn") sum+=$3; print sum}' Length_TE_Lo7_SUPERFamLevel.tsv |tail -n1
### 5376831920
gawk '{ if ($0!~"chrUn" && $0!~"unknown") sum+=$3; print sum}' Length_TE_Lo7_SUPERFamLevel.tsv |tail -n1
### 5375899302

###### TOTAL GENOME LENGTH ALL FEATURES (match_part overlapping repeat_region):
bedtools merge -i <(grep '^chr' $OUTPUT/Lo7_RefSeq_clariTE.gff3 |grep -v '^chrUn' |cut -f1,4,5 |sort -k1,1 -k2,2n) |gawk '{ sum+=$3-$2; print sum }' |tail -n1
### 5561967338 (diff avec repeat_region=225)

/home/palasser/projects/wheatomics_wp1/annotation/TE/bin/length_Family_TE.py $OUTPUT/Lo7_RefSeq_clariTE.gff3 |sort -k1,1 -k4,4n > $OUTPUT/Lo7_RefSeq_clariTE_For_Length_TE_calculation.gff3

for chr in "1R" "2R" "3R" "4R" "5R" "6R" "7R" "Un";
do
    chrom='chr'$chr
    # grep ^$chrom $OUTPUT/Lo7_RefSeq_clariTE_For_Length_TE_calculation.gff3 |gawk -v OFS='\t' '{print $NF,$5-$4}' |sed -E 's/\.[0-9]*//' \
    # |gawk -v grp=$chrom -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print grp,i,a[i]}'  |sort -k1,1 >> Length_TE_Lo7_FamLevel.tsv

    grep ^$chrom $OUTPUT/Lo7_RefSeq_clariTE_For_Length_TE_calculation.gff3 |gawk -v OFS='\t' '{print $NF,$5-$4}' |sed -E 's/_.*\t/\t/' \
    |gawk -v grp=$chrom -v OFS='\t' '{a[$1]+=$2}END{for(i in a) print grp,i,a[i]}'  |sort -k1,1 >> Length_TE_Lo7_SUPERFamLevel.tsv

    grep ^$chrom $OUTPUT/Lo7_RefSeq_clariTE.gff3 |grep '\srepeat_region\s' |cut -f9 |sed -E 's/ID=.*Family:(.*)[_|w].* Matching.*/\1/' \
    |sort |uniq -c |gawk -v grp=$chrom -v OFS='\t' '{print grp,$2,$1}' >> Nb_TE_Lo7_SUPERFamLevel.tsv
done
