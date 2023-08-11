#!/bin/bash

#SBATCH --job-name=AeU-TE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
##SBATCH --array=0-7
#SBATCH --mem=8G
#SBATCH --partition=normal

module load bedtools/2.27.1 exonerate samtools GenomeTools bioperl gdecTools

####################################################################################
##### DEFINITION VARIABLES

GENOME='/home/palasser/projects/soutien_bioinfo/UPINDER_Gill/data/FINAL_AeU_asm.fasta'
OUTPUT='/home/palasser/projects/soutien_bioinfo/UPINDER_Gill/results/annotTE_AeU'

####################################################################################
chromosomes=("1U" "2U" "3U" "4U" "5U" "6U" "7U" "Un")
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}
chrom_prefix='Chr'
chrom=${chrom_prefix}${chr}
mkdir -p $OUTPUT/${chrom}
cd $OUTPUT/${chrom}

###### Chunks
## create one fasta per chromosome
if [ $chrom = "ChrUn" ];
    then 
    gawk -v seq='ptg' 'BEGIN { RS=">" } { if ($0 ~ seq) print RS "ChrUn:"$0 }' ${GENOME} > $OUTPUT/${chrom}/${chrom}.fasta
else
    samtools faidx ${GENOME} $chrom > $OUTPUT/${chrom}/${chrom}.fasta
fi

##other methods
##gawk -v seq='chr'$chrom 'BEGIN { RS=">" } { if ($0 ~ seq) print RS $0 }' ${GENOME} > $OUTPUT/${chrom}/${chrom}.fasta
##echo $chrom > chrom.txt; seqtk subseq ${GENOME} chrom.txt > $OUTPUT/${chrom}/${chrom}.fasta

samtools faidx $OUTPUT/$chrom/${chrom}.fasta
bedtools makewindows -g $OUTPUT/${chrom}/${chrom}.fasta.fai -w 10000000 > $OUTPUT/${chrom}/${chrom}.windows.bed

bedtools getfasta -bed $OUTPUT/${chrom}/${chrom}.windows.bed -fi $OUTPUT/${chrom}/${chrom}.fasta > $OUTPUT/${chrom}/${chrom}.windows.fasta
fastaexplode -f $OUTPUT/${chrom}/${chrom}.windows.fasta -d $OUTPUT/${chrom}

files=($(find $OUTPUT/${chrom}/${chrom}*.fa))

k=$((${#files[@]}-1))   #recupere le nb d'item dans la liste 'files'

###### RepeatMasker
RMjid=`sbatch --wait --parsable -J ${chr}_RM --array=0-${k}%20 /home/palasser/projects/soutien_bioinfo/UPINDER_Gill/bin/repeatMasker_AeU.sh ${chrom}`
echo "$RMjid : RepeatMasker on ${chrom}"

if [ $(find *.fa.out.xm |wc -l) = $(($k+1)) ];
    then echo ${chrom}": all RepeatMasker output files there";
    else echo ${chrom}": missing "$(($(find $OUTPUT/${chrom}/*.fa |wc -l) - $(find $OUTPUT/${chrom}/*.fa.out.xm |wc -l) ))" RepeatMasker output files";
fi

# ## rescue example:
# # sbatch -J 6A_RMrescue -p normal -c 16 --mem=16G --wrap="module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker; cd /home/palasser/projects/soutien_bioinfo/UPINDER_Gill/results/annotTE_AeU/Chr1U; RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa 16 -q /home/palasser/projects/soutien_bioinfo/UPINDER_Gill/results/annotTE_AeU/Chr1U/Chr1U:80000000-90000000.fa"


###### clariTE
## option `--wait` pour que le script principal attende la fin du job array clariTE avant cmd suivante
CLARITEjid=`sbatch --wait --parsable -J ${chr}_clariTE --dependency=afterok:$RMjid --array=0-$k%20 /home/palasser/projects/soutien_bioinfo/UPINDER_Gill/bin/clariTE_AeU.sh ${chrom}`

echo "$CLARITEjid : clariTE on ${chrom}"

if [ $(find *.fa.out_anno.embl |wc -l) = $(($k+1)) ];
    then echo ${chrom}": all clariTE output files there";
    else echo ${chrom}": missing "$(($(find $OUTPUT/${chrom}/*.fa |wc -l) - $(find $OUTPUT/${chrom}/*.fa.out_anno.embl |wc -l) ))" clariTE output files";
fi

### /!\ creation d'un script specifique : /storage/groups/gdec/bin/scripts/embl2gff.AeU.pl >>> PAS PROPRE A REVOIR
###### embl2gff: takes all chunks embl and produces one merged .gff
## le script embl2gff.AeU.pl numerote les TE par chrom (et non par chunk) en repartant de zero pour chaque chrom a condition
## de fournir la liste ordonnee des fichiers pour l'ensemble des chunks par chromosome
embl_files=$(\ls -1 $OUTPUT/${chrom}/*.embl |sort -t ':' -k2,2n |tr -s '\n' ' ')
/storage/groups/gdec/bin/scripts/embl2gff.AeU.pl -renan -note -l 10 $embl_files > $OUTPUT/AeU_RefSeq_${chr}_clariTE.gff

###### gff to gff3: recalcul des coordonnees avec gawk (transfo de relatives aux chunks vers relatives au chrom) + mise en forme gff3 par GenomeTools sort
grep -v $'\t''region' $OUTPUT/AeU_RefSeq_${chr}_clariTE.gff |sed 's/embl2gff.AeU.pl-2.1/clariTE/' \
|gawk -v FS='\t' -v OFS='\t' '/^Chr/ {match($1, /:[0-9]+-/); $4=$4+substr($1, RSTART+1, RLENGTH-2); $5=$5+substr($1, RSTART+1, RLENGTH-2); print}' \
|gt gff3 -sort -tidy -retainids 1> $OUTPUT/AeU_RefSeq_${chr}_clariTE.gff3  2> $OUTPUT/AeU_RefSeq_${chr}_clariTE_gt.log

###### friendly gff3
if [ $chrom = "ChrUn" ];
    then 
    endchrom=$(gawk '{sum+=$2; print sum}' $OUTPUT/${chrom}/${chrom}.fasta.fai |tail -n1)
    #le premier grep garde la premiere "##sequence-region" et les annotation de TEs
    grep -P '##sequence-region *'$chrom':ptg000043l|^Chr|###' $OUTPUT/AeU_RefSeq_${chr}_clariTE.gff3 \
    |gawk -v LG=$chrom -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"LG"\t1\t"end; print }' \
    |sed -E 's/'$chrom':[0-9]*-[0-9]*/'$chrom'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > $OUTPUT/AeU_RefSeq_${chr}_clariTE_friendly.gff3
else
    endchrom=$(grep $chrom ${GENOME}.fai |cut -f2)
    #le premier grep retire toute les "##sequence-region" sauf la premiere
    grep -v -P '##sequence-region *'$chrom':[1-9]' $OUTPUT/AeU_RefSeq_${chr}_clariTE.gff3 \
    |gawk -v LG=$chrom -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"LG"\t1\t"end; print }' \
    |sed -E 's/'$chrom':[0-9]*-[0-9]*/'$chrom'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > $OUTPUT/AeU_RefSeq_${chr}_clariTE_friendly.gff3
fi


echo "Check clariTE_gff3 validity for ${chrom}"
gt gff3validator $OUTPUT/AeU_RefSeq_${chr}_clariTE_friendly.gff3


rm $OUTPUT/${chrom}/${chrom}.fasta $OUTPUT/${chrom}/${chrom}.fasta.fai $OUTPUT/${chrom}/${chrom}.windows.fasta $OUTPUT/${chrom}/${chrom}.windows.bed