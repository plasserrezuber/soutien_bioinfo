
#!/bin/bash
#SBATCH --job-name=mulifa_split
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1
#SBATCH --array=1-94

#################################################################################
OUTPUT='/home/newatil/multifa_split'
mkdir -p ${OUTPUT}
### errors management
### set -x #mode debug on
set -o errexit #pour que le script s'arrete en cas d erreur ignoree
set -o nounset #force l initialisation des variables
IFS=$'\n\t'

### modules management
module load gcc/8.1.0 seqtk/1.3 

########################################################################################################
fifa=$(\ls -1 /home/newatil/pbaa/*passed_cluster_sequences.fasta|sed -n ${SLURM_ARRAY_TASK_ID}p)
fifai=$(\ls -1 /home/newatil/pbaa/*passed_cluster_sequences.fasta.fai|sed -n ${SLURM_ARRAY_TASK_ID}p)

echo '*****************************************************************************'
echo "start job  ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} on ${HOSTNAME} at `date`"

######## seqtk split multifasta on single fasta #############################################################
echo '######## seqtk runing ################'
for seq in $(cut -f1 $fifai); do seqtk subseq $fifa <(echo $seq) > ${OUTPUT}/$seq.fasta; done

echo '################# DONE ########################'
