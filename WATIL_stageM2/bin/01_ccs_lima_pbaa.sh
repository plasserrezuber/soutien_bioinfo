#!/bin/bash
#SBATCH --job-name=pbaa
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p normal
#SBATCH --cpus-per-task=8
#SBATCH --array=1-94

################################################################################
### errors management
### set -x #mode debug on
set -o errexit #pour que le script s'arrete en cas d erreur ignoree
set -o nounset #force l initialisation des variables
IFS=$'\n\t'

### modules management
module load gcc/8.1.0 smrttools/10.1.0.119588 samtools/1.9


########################################################################################################

SCRATCHDIR=/storage/scratch/$USER/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p -m 750 ${SCRATCHDIR}
cd ${SCRATCHDIR}
OUTPUT='/home/newatil/pbaa'
mkdir -p ${OUTPUT}

########################################################################################################
### partie realisee par Vincent Pailler
#ccs -j 128 --min-passes 3 *.subreads.bam ccs.bam

### partie realisee par Nezha Watil
#lima --split-bam-named -j 32 --ccs --peek-guess  /home/newatil/data/ccs.bam /home/newatil/data/SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta /home/newatil/PacBio/demux.bam
########################################################################################################

FILE=$(\ls -1 /home/newatil/PacBio/demux.*.bam |sed -n ${SLURM_ARRAY_TASK_ID}p)
SHORTNAME=$(basename ${FILE} |cut -d'.' -f1-2)


echo '*****************************************************************************'
echo "start job  ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} on ${HOSTNAME} at `date`"
echo '************ PacBio Amplicon Analysis of sample' ${SHORTNAME}' ***************'

######## smrttools/bam2fastq #############################################################
echo '######## convert bam to fastq ################'
samtools fastq ${FILE} > ${SHORTNAME}.fastq

echo '######## create indedx with fqidx #############'

samtools fqidx ${SHORTNAME}.fastq

echo '################ Runing pbaa ##################################'

pbaa cluster -j 8 --max-consensus-reads 1000 --log-file ${SHORTNAME}.log /home/palasser/soutien_bioinfo/WATIL_stageM2/18_amplicons_GLU_CSRefSeqv1.fasta ${SHORTNAME}.fastq pbaa_${SHORTNAME}

############ move output data ###################
echo '####################### start moving data #####################################'

rm ${SCRATCHDIR}/${SHORTNAME}.fastq
rm ${SCRATCHDIR}/${SHORTNAME}.fastq.fai
mv ${SCRATCHDIR}/${SHORTNAME}.log ${OUTPUT}/
mv ${SCRATCHDIR}/pbaa_* ${OUTPUT}/ && rm -Rf ${SCRATCHDIR}

cd ${OUTPUT}
sed -i 's/HMW_/HMW/g' pbaa_${SHORTNAME}_passed_cluster_sequences.fasta
sed -i 's/HMW_/HMW/g' pbaa_${SHORTNAME}_passed_cluster_sequences.fasta.fai

echo '********* DONE *******************'

