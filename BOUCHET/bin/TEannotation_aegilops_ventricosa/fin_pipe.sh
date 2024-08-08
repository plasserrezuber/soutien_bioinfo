#!/bin/bash

#SBATCH --job-name=embl2gff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --array=0-373%20
#SBATCH --mem=8G
#SBATCH --partition=fast
#SBATCH --qos=fast

module load GenomeTools/1.5.5 bioperl/1.7.0_rc5

OUTPUT='/home/palasser/projects/soutien_bioinfo/BOUCHET/results/annotTE_aegilops_ventricosa/chrom/ctg1Mb'
cd $OUTPUT

### NB: pas de TE dans contig ptg000372l
contigs=($(cut -f1 ../ctg1Mb.fasta.fai |grep -v 'ptg000372l'))

c=${contigs[$SLURM_ARRAY_TASK_ID]}

files=$(find ${c}*embl |sort -t ':' -k2,2n |tr -s '\n' ' ')
/home/palasser/projects/soutien_bioinfo/BOUCHET/bin/TEannotation_aegilops_ventricosa/bin/embl2gff.pl -RMclariTE -featurePrefix ae_ventricosa_ -source clariTE -note -l 10 ${files} > ${c}_clariTE.gff


###### gff to gff3: new coordinates calculation with gawk commands "match" + "substr" (from chunks relative to chrom relative)
grep -v $'\t''region' ${c}_clariTE.gff \
|awk -v FS='\t' -v OFS='\t' '/^ptg/ {match($1, /:[0-9]+-/); $4=$4+substr($0, RSTART+1, RLENGTH-2); $5=$5+substr($0, RSTART+1, RLENGTH-2); print}' \
|gt gff3 -sort -tidy -retainids 1> ${c}_clariTE_tmp.gff3

### variable $endchrom to give the right sequence-region coordinate in gff3
endchrom=$(grep ${c} ../ctg1Mb.fasta.fai |cut -f2)

# grep -v command keeps only one "sequence-region" line starting with zero coordinate
# gawk and sed commands to format in a more friendly way
grep -v -P '##sequence-region *'$c':[1-9]' ${c}_clariTE_tmp.gff3 \
|awk -v LG=$c -v end=$endchrom 'BEGIN{FS="\t";OFS="\t"} { if ($0~"##sequence-region") $0="##sequence-region\t"LG"\t1\t"end; print }' \
|sed -E 's/'$c':[0-9]*-[0-9]*/'$c'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > ${c}_clariTE.gff3
