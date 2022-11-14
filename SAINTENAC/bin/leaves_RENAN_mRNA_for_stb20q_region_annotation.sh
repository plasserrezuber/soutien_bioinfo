#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1

ml samtools

support_dir="/home/palasser/projects/soutien_bioinfo/SAINTENAC/RNAseq_leaves_Renan_Wheatomics"
mkdir -p ${support_dir}
dir="/home/palasser/projects/wheatomics_wp1/annotation/transcripts/results/"

for name in CGC_AAHN CGC_AAHO CGC_AAIH CGC_AAII CGC_AAJB CGC_AAJC;
do
    basename_gtf=$(basename $(find ${dir}/${name}*_q60.gtf))
    basename_bam=$(basename $(find ${dir}/${name}*_sorted_q60.bam))
    ln -s ${dir}/${name}*_q60.gtf ${support_dir}/${basename_gtf}
    ln -s ${dir}/${name}*_sorted_q60.bam ${support_dir}/${basename_bam}
    #samtools index -c ${support_dir}/${name}*_sorted_q60.bam
done