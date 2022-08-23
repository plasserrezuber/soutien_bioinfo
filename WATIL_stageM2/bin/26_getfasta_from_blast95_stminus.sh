#!/bin/bash
#SBATCH --job-name=getfasta
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH -p debug
##SBATCH --qos=fast
#SBATCH --cpus-per-task=4
#SBATCH --array=1-1404


#################################################################################
OUTPUT='/home/newatil/fasta_files_blast95'
mkdir -p ${OUTPUT}

### modules management
module load seqtk/1.3 gdecTools bioperl 


########################################################################################################

FILE=$(\ls -1 /home/newatil/blastn_pid95/demux.bc*test.fasta_stminus.blastn |sed -n ${SLURM_ARRAY_TASK_ID}p)
SHORTNAME=$(basename ${FILE} |cut -d'.' -f1-2)


echo '*****************************************************************************'

##########################################################################################################

for f in ${FILE};
do
  seqid95=$(awk '{ print $2}' ${f});
  for s in ${seqid95} ;
  do
    fgrep ${s} /home/newatil/testPY/${SHORTNAME}.fasta|cut -f 1 -d ' ' | sed 's/>//g' >> ${SHORTNAME}_stminus.list
    seqtk subseq /home/newatil/testPY/${SHORTNAME}.fasta ${SHORTNAME}_stminus.list > ${SHORTNAME}_stminus.fasta
    revcom.pl -o ${OUTPUT} ${SHORTNAME}_stminus.fasta
  done
done

echo '************************************DONE***************************************'
