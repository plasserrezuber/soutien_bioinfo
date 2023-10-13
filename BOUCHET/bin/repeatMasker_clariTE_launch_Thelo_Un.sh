#!/bin/bash

#SBATCH --job-name=annoTE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --partition=gdec

module load bedtools/2.27.1 exonerate/2.4.0 samtools/1.3 GenomeTools/1.5.5 bioperl/1.7.0_rc5 gdecTools/1.1

####################################################################################
##### DEFINITION VARIABLES

GENOME=/storage/groups/gdec/shared/thinopyrum_elongatum/GCA_011799875.1_ASM1179987v1/GCA_011799875.1_ASM1179987v1_genomic.fna
OUTPUT=/home/palasser/projects/soutien_bioinfo/BOUCHET/results/annotTE_Thelo

####################################################################################
chr='Un'
chrom='chr'$chr
mkdir -p $OUTPUT/${chrom}
cd $OUTPUT/${chrom}

# ###### Chunks
## create one fasta per chromosome
if [ $chrom = "chrUn" ];
    then 
    gawk -v seq='JAAAX' 'BEGIN { RS=">" } { if ($0 ~ seq) print RS "chrUn:"$0 }' ${GENOME} > ${chrom}.fasta
else
    samtools faidx $GENOME $accessNum |sed 's/>'$accessNum'/>'$chrom'/' > ${chrom}.fasta
fi

samtools faidx ${chrom}.fasta
# bedtools makewindows -g ${chrom}.fasta.fai -w 10000000 > ${chrom}.windows.bed

# bedtools getfasta -bed ${chrom}.windows.bed -fi ${chrom}.fasta > ${chrom}.windows.fasta
# fastaexplode -f ${chrom}.windows.fasta -d $OUTPUT/${chrom}

files=($(find ${chrom}*.fa))
k=$((${#files[@]}-1))  #recupere le nb d'item dans la liste 'files'


# ###### RepeatMasker
# RMjid=`sbatch --wait --parsable -J ${chr}_RM --array=0-${k}%20 /home/palasser/projects/soutien_bioinfo/BOUCHET/bin/repeatMasker_Thelo.sh ${chrom}`

# echo "$RMjid : RepeatMasker on ${chrom}"

# if [ $(find ${chrom}:*.fa.out.xm |wc -l) = $(($k+1)) ];
#     then echo ${chrom}": all RepeatMasker output files there";
#     else echo ${chrom}": missing "$(($(find ${chrom}:*.fa |wc -l) - $(find ${chrom}:*.fa.out.xm |wc -l) ))" RepeatMasker output files";
# fi

# ## rescue example:
# # sbatch -J 6A_RMrescue -p normal -c 16 --mem=8G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/projects/soutien_bioinfo/BOUCHET/results/annotTE_Thelo/chr6A; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/projects/soutien_bioinfo/BOUCHET/results/annotTE_Thelo/chr6A/chr6A:80000000-90000000.fa"
# # sbatch -J 7B_RMrescue -p normal -c 16 --mem=8G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/projects/soutien_bioinfo/BOUCHET/results/annotTE_Thelo/chr7B; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/projects/soutien_bioinfo/BOUCHET/results/annotTE_Thelo/chr7B/chr7B:440000000-450000000.fa"


###### clariTE
## option `--wait` pour que le script principal attende la fin du job array clariTE avant cmd suivante
# --dependency=afterok:$RMjid 
# CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --array=0-$k%20 /home/palasser/projects/soutien_bioinfo/BOUCHET/bin/clariTE_Thelo.sh ${chrom}`
# # CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --array=0-$k%20 /home/palasser/projects/soutien_bioinfo/BOUCHET/bin/clariTE_Thelo.sh ${chrom}`

# echo "$CLARITEjid : clariTE on ${chrom}"

# if [ $(find ${chrom}:*.fa.out_anno.embl |wc -l) = $(($k+1)) ];
#     then echo ${chrom}": all clariTE output files there";
#     else echo ${chrom}": missing "$(($(find ${chrom}:*.fa |wc -l) - $(find ${chrom}:*.fa.out_anno.embl |wc -l) ))" clariTE output files";
# fi

###### embl2gff: takes all chunks embl and produces one merged .gff
## le script embl2gff.generic.pl numerote les TE par chrom (et non par chunk) en repartant de zero pour chaque chrom a condition
## de fournir la liste ordonnee des fichiers pour l'ensemble des chunks par chromosome
embl_files=$(\ls -1 ${chrom}*.embl |sort -t ':' -k2,2n |tr -s '\n' ' ')
/storage/groups/gdec/bin/scripts/embl2gff.pl -RMclariTE -featurePrefix Thelo_chrUn -note -l 10 $embl_files > $OUTPUT/Thelo_RefSeq_${chr}_clariTE.gff

###### gff to gff3: recalcul des coordonnees avec gawk (transfo de relatives aux chunks vers relatives au chrom) + mise en forme gff3 par GenomeTools sort
grep -v $'\t''region' $OUTPUT/Thelo_RefSeq_${chr}_clariTE.gff |sed 's/embl2gff.pl-2.2/clariTE/' \
|gawk -v FS='\t' -v OFS='\t' '/^chr/ {match($1, /:[0-9]+-/); $4=$4+substr($0, RSTART+1, RLENGTH-2); $5=$5+substr($0, RSTART+1, RLENGTH-2); print}' \
|gt gff3 -sort -tidy -retainids 1> $OUTPUT/Thelo_RefSeq_${chr}_clariTE.gff3  2> $OUTPUT/Thelo_RefSeq_${chr}_clariTE_gt.log


###### friendly gff3
if [ $chrom = "chrUn" ];
    then 
    endchrom=$(gawk '{sum+=$2; print sum}' $OUTPUT/${chrom}/${chrom}.fasta.fai |tail -n1)
    #le premier grep garde la premiere "##sequence-region" et les annotation de TEs
    grep -P '##sequence-region *'$chrom':JAAAXO010000008|^chr|###' $OUTPUT/Thelo_RefSeq_${chr}_clariTE.gff3 \
    |gawk -v LG=$chrom -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"LG"\t1\t"end; print }' \
    |sed -E 's/'$chrom':JAAAXO.*embl.1/'$chrom'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > $OUTPUT/Thelo_RefSeq_${chr}_clariTE_friendly.gff3
else
    endchrom=$(grep $accessNum ${GENOME}.fai |cut -f2)
    grep -v -P '##sequence-region *'$chrom':[1-9]' $OUTPUT/Thelo_RefSeq_${chr}_clariTE.gff3 \
    |gawk -v chr=$chrom -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"chr"\t1\t"end; print }' \
    |sed -E 's/'$chrom':[0-9]*-[0-9]*/'$chrom'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > $OUTPUT/Thelo_RefSeq_${chr}_clariTE_friendly.gff3
fi

echo "Check clariTE_gff3 validity for ${chrom}"
gt gff3validator $OUTPUT/Thelo_RefSeq_${chr}_clariTE_friendly.gff3


rm ${chrom}.fasta ${chrom}.fasta.fai ${chrom}.windows.fasta ${chrom}.windows.bed

