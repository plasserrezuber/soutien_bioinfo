#!/bin/bash
#SBATCH --job-name=cd-hit
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --array=1-1694


#################################################################################
OUTPUT='/home/newatil/cdhit_mergedfasta'

### errors management
### set -x #mode debug on
set -o errexit #pour que le script s'arrete en cas d erreur ignoree
set -o nounset #force l initialisation des variables
IFS=$'\n\t'

### modules management
module load cdhit/4.6.8

### working directory analysis in scratch
SCRATCHDIR=/storage/scratch/$USER/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p -m 750 ${SCRATCHDIR}
cd ${SCRATCHDIR}
mkdir -p ${OUTPUT}

########################################################################################################

FILE=$(\ls -1 /home/newatil/fasta_files_merged/demux*merged.fasta|sed -n ${SLURM_ARRAY_TASK_ID}p)
SHORTNAME=$(basename ${FILE} |cut -d'.' -f1-2)


echo '*****************************************************************************'
echo "start job  ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} on ${HOSTNAME} at `date`"
echo '************ PacBio Amplicon Analysis of sample' ${SHORTNAME}' ***************'

########   Runing cd-hit program #############################################################
echo '######## convert bam to fastq ################'
cd-hit-est -T 8 -M 16000 -g 1 -c 0.99 -aL 0.9 -aS 0.5 -d 0 -i ${FILE} -o ${SHORTNAME}.fasta

echo '#######################################################'

mv ${SCRATCHDIR}/* ${OUTPUT}/ && rm -Rf ${SCRATCHDIR}

echo '********* DONE *******************'
