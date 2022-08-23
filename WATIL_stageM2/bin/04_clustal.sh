#!/bin/bash
#SBATCH --job-name=clustalo
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --array=0-17
#################################################################################
echo '******************************************'
echo 'Alignement multiple'
echo '******************************************'

### modules management
module load ClustalOmega/1.2.4

### working directory analysis

OUTPUT='/home/newatil/pbaa_rc_clustalo'
mkdir -p ${OUTPUT}
cd ${OUTPUT}

########################################################################################################
genes=("HMW1Ax" "HMW1Ay" "HMW1Bx" "HMW1By" "HMW1Dx" "HMW1Dy" "LMW1A_i_1" "LMW1A_i_2" "LMW1A_m_3" "LMW1B_m_4" "LMW1B_m_5" "LMW1D_m_1" "LMW1D_m_3" "LMW1D_m_4" "LMW1D_m_5" "LMW1D_m_6" "LMW1D_m_7" "LMW1D_m_8")

g=${genes[$SLURM_ARRAY_TASK_ID]}
cat /home/newatil/blastn_pbaa_rc/sample*${g}*fasta > pbaa_${g}.fasta

if [ $g = "LMW1D_m_7" ];
	then
	gawk 'BEGIN {RS=">" ; FS="\n"; ORS=">"} {if ($1!~"bc2004" && $1!~"bc2010" && $1!~"bc2030") print}' pbaa_${g}.fasta |head -n -1 \
	> pbaa_${g}_curated.fasta

elif [ $g = "LMW1D_m_4" ];
  then
	gawk 'BEGIN {RS=">" ; FS="\n"; ORS=">"} {if ($1!~"bc2048") print}' pbaa_${g}.fasta |head -n -1 \
  > pbaa_${g}_curated.fasta

elif [ $g = "LMW1D_m_5" ];
  then
  gawk 'BEGIN {RS=">" ; FS="\n"; ORS=">"} {if ($1!~"bc2004" && $1!~"bc2050" && $1!~"bc2048") print}' pbaa_${g}.fasta |head -n -1 \
  > pbaa_${g}_curated.fasta

elif [ $g = "LMW1B_m_5" ];
  then
  gawk 'BEGIN {RS=">" ; FS="\n"; ORS=">"} {if ($1!~"bc2023" && $1!~"bc2069" && $1!~"bc2004" && $1!~"bc2005") print}' pbaa_${g}.fasta |head -n -1 \
  > pbaa_${g}_curated.fasta

else
	mv pbaa_${g}.fasta pbaa_${g}_curated.fasta
fi

######## Run the program

clustalo --threads 8 -i pbaa_${g}_curated.fasta -t DNA -o pbaa_${g}_aling.fasta
###########################################################################################################
echo '********************************DONE********************************'
