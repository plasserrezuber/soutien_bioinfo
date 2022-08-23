#!/bin/bash
#SBATCH --job-name=emboss_cons
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -p debug
#SBATCH --cpus-per-task=2
#SBATCH --array=1-1279

OUTPUT='/home/newatil/Cons_for_merged_files'
mkdir -p ${OUTPUT}
cd ${OUTPUT}

ml java/oracle-1.8.0_45 emboss/6.5.7

FILE=$(\ls -1 /home/newatil/Cons_for_merged_files/demux.bc20*test_merged.clustal|sed -n ${SLURM_ARRAY_TASK_ID}p)
SHORTNAME=$( basename ${FILE}|cut -d '_' -f1-3)

cons -sequence ${FILE} -outseq ${SHORTNAME}.cons.fasta -plurality 1 -setcase 0 -name ${SHORTNAME}

echo '************************* finish *********************************'
