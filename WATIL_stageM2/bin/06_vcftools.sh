#!/bin/bash
#SBATCH --job-name=mulifa_split
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1
#SBATCH --array=1-18

ml java/oracle-1.8.0_45 vcftools/v0.1.16

## fichiers vcf obtenus en local avec msa2vcf copies dans /home/newatil/msa2vcf/
############################################################################
################## statistical analysis of glutinins diversity #############
############################################################################

f=$(\ls -1 /home/newatil/msa2vcf/pbaa*aling.fasta.vcf |sed -n ${SLURM_ARRAY_TASK_ID}p)

#################  allelic frequencies #####################################

vcftools --vcf ${f} --freq  --out ${f} # renvoie la frequences allélique en %
vcftools --vcf ${f} --counts  --out ${f} #renvoie le nombre des alléles

############### Tajima D ###################################################

vcftools --vcf ${f} --TajimaD 1000 --out ${f}

############### indels STATISTICS #############################################

vcftools --vcf ${f} --hist-indel-len --out ${f} #renvoie nombre des SNP,insertions,deletions

############### HWE STATISTICS #############################################

vcftools --vcf ${f} --hardy  --out ${f} #Reports a p-value for each site from a Hardy-Weinberg Equilibrium test

################# heterozygosity rate ###################################### 

vcftools --vcf ${f} --het  --out ${f} #Calculates a measure of heterozygosity

################# TRANSITION/TRANSVERSION STATISTICS #######################

vcftools --vcf ${f} --TsTv-summary  --out ${f} #Calculates a simple summary of all Transitions and Transversions
