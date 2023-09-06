#!/bin/bash

#SBATCH --job-name=annoTE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --array=0-7
#SBATCH --mem=8G
#SBATCH --partition=normal

module load bedtools/2.27.1 exonerate/2.4.0 samtools/1.3 GenomeTools/1.5.5 bioperl/1.7.0_rc5 gdecTools/1.1

####################################################################################
##### DEFINITION VARIABLES

DATA='/storage/groups/gdec/shared/Secale_cereale/Lo7_RefSeq'    #Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta
OUTPUT='/home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7'

####################################################################################
chromosomes=("1R" "2R" "3R" "4R" "5R" "6R" "7R" "Un")
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}
chrom='chr'$chr
mkdir -p $OUTPUT/${chrom}
cd $OUTPUT/${chrom}

###### Chunks
## create one fasta per chromosome
samtools faidx $DATA/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta $chrom > $OUTPUT/${chrom}/${chrom}.fasta
##other methods
##gawk -v seq='chr'$chrom 'BEGIN { RS=">" } { if ($0 ~ seq) print RS $0 }' $DATA/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta > $OUTPUT/${chrom}/${chrom}.fasta
##echo $chrom > chrom.txt; seqtk subseq $DATA/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta chrom.txt > $OUTPUT/${chrom}/${chrom}.fasta

samtools faidx $OUTPUT/$chrom/${chrom}.fasta
bedtools makewindows -g $OUTPUT/${chrom}/${chrom}.fasta.fai -w 10000000 > $OUTPUT/${chrom}/${chrom}.windows.bed

bedtools getfasta -bed $OUTPUT/${chrom}/${chrom}.windows.bed -fi $OUTPUT/${chrom}/${chrom}.fasta > $OUTPUT/${chrom}/${chrom}.windows.fasta
fastaexplode -f $OUTPUT/${chrom}/${chrom}.windows.fasta -d $OUTPUT/${chrom}

files=($(find $OUTPUT/${chrom}/${chrom}*.fa))
k=$((${#files[@]}-1))  #recupere le nb d'item dans la liste 'files'


###### RepeatMasker
RMjid=`sbatch --wait --parsable -J ${chr}_RM --array=0-${k}%20 /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/bin/repeatMasker_Lo7.sh ${chrom}`

echo "$RMjid : RepeatMasker on ${chrom}"

if [ $(find ${chrom}:*.fa.out.xm |wc -l) = $(($k+1)) ];
    then echo ${chrom}": all RepeatMasker output files there";
    else echo ${chrom}": missing "$(($(find $OUTPUT/${chrom}/${chrom}:*.fa |wc -l) - $(find $OUTPUT/${chrom}/${chrom}:*.fa.out.xm |wc -l) ))" RepeatMasker output files";
fi

# ## rescue example:
# # sbatch -J 6A_RMrescue -p normal -c 16 --mem=8G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7/chr6A; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7/chr6A/chr6A:80000000-90000000.fa"
# # sbatch -J 7B_RMrescue -p normal -c 16 --mem=8G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7/chr7B; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7/chr7B/chr7B:440000000-450000000.fa"


###### clariTE
## option `--wait` pour que le script principal attende la fin du job array clariTE avant cmd suivante
CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --dependency=afterok:$RMjid --array=0-$k%20 /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/bin/clariTE_Lo7.sh ${chrom}`
# CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --array=0-$k%20 /home/palasser/projects/soutien_bioinfo/BOUGUENNEC/bin/clariTE_Lo7.sh ${chrom}`

echo "$CLARITEjid : clariTE on ${chrom}"

if [ $(find ${chrom}:*.fa.out_anno.embl |wc -l) = $(($k+1)) ];
    then echo ${chrom}": all clariTE output files there";
    else echo ${chrom}": missing "$(($(find $OUTPUT/${chrom}/${chrom}:*.fa |wc -l) - $(find $OUTPUT/${chrom}/${chrom}:*.fa.out_anno.embl |wc -l) ))" clariTE output files";
fi

###### embl2gff: takes all chunks embl and produces one merged .gff
## le script embl2gff.generic.pl numerote les TE par chrom (et non par chunk) en repartant de zero pour chaque chrom a condition
## de fournir la liste ordonnee des fichiers pour l'ensemble des chunks par chromosome
embl_files=$(\ls -1 $OUTPUT/${chrom}/${chrom}*.embl |sort -t ':' -k2,2n |tr -s '\n' ' ')
/storage/groups/gdec/bin/scripts/embl2gff.generic.pl -renan -note -l 10 $embl_files > $OUTPUT/Lo7_RefSeq_${chr}_clariTE.gff

###### gff to gff3: recalcul des coordonnees avec gawk (transfo de relatives aux chunks vers relatives au chrom) + mise en forme gff3 par GenomeTools sort
grep -v $'\t''region' $OUTPUT/Lo7_RefSeq_${chr}_clariTE.gff |sed 's/embl2gff.generic.pl-2.1/clariTE/' \
|gawk -v FS='\t' -v OFS='\t' '/^chr/ {match($1, /:[0-9]+-/); $4=$4+substr($0, RSTART+1, RLENGTH-2); $5=$5+substr($0, RSTART+1, RLENGTH-2); print}' \
|gt gff3 -sort -tidy -retainids 1> $OUTPUT/Lo7_RefSeq_${chr}_clariTE.gff3  2> $OUTPUT/Lo7_RefSeq_${chr}_clariTE_gt.log

###### friendly gff3
endchrom=$(grep $chrom $DATA/Secale_cereale_Lo7_RefSeq_RabanusWallace_etal_2020_pseudomolecules.fasta.fai |cut -f2)
grep -v -P '##sequence-region *'$chrom':[1-9]' $OUTPUT/Lo7_RefSeq_${chr}_clariTE.gff3 \
|gawk -v chr=$chrom -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"chr"\t1\t"end; print }' \
|sed -E 's/'$chrom':[0-9]*-[0-9]*/'$chrom'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > $OUTPUT/Lo7_RefSeq_${chr}_clariTE_friendly.gff3

sed -i 's/=TraesRe_chr/=SecerLo7_chr/g' Lo7_RefSeq_${chr}_clariTE_friendly.gff3

echo "Check clariTE_gff3 validity for ${chrom}"
gt gff3validator $OUTPUT/Lo7_RefSeq_${chr}_clariTE_friendly.gff3


rm $OUTPUT/${chrom}/${chrom}.fasta $OUTPUT/${chrom}/${chrom}.fasta.fai $OUTPUT/${chrom}/${chrom}.windows.fasta $OUTPUT/${chrom}/${chrom}.windows.bed

