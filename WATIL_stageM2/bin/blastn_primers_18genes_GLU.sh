#!/bin/bash

ml ncbi-blast/2.11.0+ gdecTools/1.1 triannotTools/1.2 bioperl/1.7.0_rc5

sbatch -p debug -c 8 --wrap="ml ncbi-blast/2.11.0+;
blastn -num_threads 8 -outfmt 6 -word_size 5 -dust no \
-query /home/palasser/soutien_bioinfo/WATIL_stageM2/primer_18_amplicons_CS_v1.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-out /home/palasser/soutien_bioinfo/WATIL_stageM2/blastn_primer_vs_CSRefSeqv1.blastn"

sed 's/TE_LMW1d_m_6_R6D/TE_LMW1D_m_6_R6D/' blastn_primer_vs_CSRefSeqv1.blastn |sed 's/TE_LMW1D-m-4-/TE_LMW1D_m_4_/' |fgrep -wv '4343974' > blastn_primer_vs_CSRefSeqv1_modif.blastn

for chr in "1A" "1B" "1D";
do
    parseBlast.pl -f -if m8  blastn_primer_vs_CSRefSeqv1_modif.blastn |egrep -w "chr"$chr |gawk -v chrom=$chr '{ if ($2~chrom && $3~chrom) print $0}'
done |sort -k1,1n > /home/palasser/soutien_bioinfo/WATIL_stageM2/blastn_primer_vs_CSRefSeqv1_parsed.blastn

gawk -v OFS='\t' '{ print $3,$11,$12,$2}' /home/palasser/soutien_bioinfo/WATIL_stageM2/blastn_primer_vs_CSRefSeqv1_parsed.blastn |sort -k1,1 -k2,2n -k3,3n \
|gawk -v OFS='\t' '{ if ($2>$3) print $1,$3,$2,$4,$4; else print $0,$4}' \
|gawk -v OFS='\t' '{$5=substr($5, 4,9); print }' |cut -d'_' -f1-7 |sed -E 's/\tHMW(.+)_[R|F].*$/\tHMW\1/' \
|awk -v OFS="\t" 'NR!=1 {start_before=a[2]; name_before=a[5]} {split($0,a,FS)} { print $1,start_before,$3,$5,name_before }' \
|gawk -v OFS='\t' '{ if ($4==$5 && NF==5) print $1,$2,$3,$4 }' \
> /home/palasser/soutien_bioinfo/WATIL_stageM2/blastn_primer_vs_CSRefSeqv1_parsed.bed

ml bedtools

bedtools getfasta -name -fi /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-bed /home/palasser/soutien_bioinfo/WATIL_stageM2/blastn_primer_vs_CSRefSeqv1_parsed.bed \
-fo /home/palasser/soutien_bioinfo/WATIL_stageM2/18_amplicons_GLU_CSRefSeqv1.fasta
