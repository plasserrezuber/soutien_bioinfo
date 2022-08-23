#!/bin/bash
#SBATCH --job-name=GetCons
#SBATCH -p normal
##BATCH --qos=fast
#SBATCH --mem=32G
#BATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --array=1-1504

OUTPUT='/home/newatil/ExtractCons'

SCRATCHDIR=/storage/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p -m 700 ${SCRATCHDIR}
cd ${SCRATCHDIR}
mkdir -p ${OUTPUT}

#########################################################################
#GENES=$(fgrep '>' /home/newatil/PacBio/E779-Primers.txt |sed 's/>Glu_//' |cut -d'-' -f1 |sort |uniq)
#for g in ${GENES},do
#FILE=$(\ls -1 /home/newatil/cd-hit_filtredseq/demux*_${g}_test.fasta.clstr|sed -n ${SLURM_ARRAY_TASK_ID}p)
#DATA=$(\ls -1 /home/newatil/blast_filtred_id90/demux*_${g}_test.fasta|sed -n ${SLURM_ARRAY_TASK_ID}p)
########################################################################



FILE=$(\ls -1 /home/newatil/cd-hit_filtredseq/demux*_test.fasta.clstr|sed -n ${SLURM_ARRAY_TASK_ID}p)
DATA=$(\ls -1 /home/newatil/blast_filtred_id90/demux*test.fasta|sed -n ${SLURM_ARRAY_TASK_ID}p)
SHORTNAME=$(basename ${FILE} |cut -d'.' -f1-2)

echo '******************************************'
echo "load required modules at `date`"
echo '******************************************'
 
module load gcc/8.1.0 python/3.7.1 seqtk/1.3 java/oracle-1.8.0_45 emboss/6.5.7 ClustalOmega/1.2.4


echo '******************************************'
echo "start runing at `date`"
echo '******************************************'


#for f in $(FILE);
#do

echo '******************************************'
echo 'Extract sequences'
echo '******************************************'


python /home/newatil/bin/extract_representativeSeq.py -i ${FILE}

seqtk subseq ${DATA} /home/newatil/cd-hit_filtredseq/${SHORTNAME}.fasta.txt > ${SHORTNAME}.seqtk.fasta
 


echo '******************************************'
echo 'Alignement multiple'
echo '******************************************'


clustalo -i ${SHORTNAME}.seqtk.fasta -t  DNA -o ${SHORTNAME}.clustal.fasta --threads 8



echo '******************************************'
echo 'Generation de consensus'
echo '******************************************'


cons -sequence  ${SHORTNAME}.clustal.fasta  -outseq ${SHORTNAME}.cons.fasta


echo '******************************************'
echo 'remove scratch'
echo '******************************************'


#mv ${SCRATCHDIR}/${SHORTNAME}.fasta.txt ${OUTPUT}/
#mv ${SCRATCHDIR}/${SHORTNAME}.seqtk.fasta ${OUTPUT}/
#mv ${SCRATCHDIR}/${SHORTNAME}.clustal.fasta ${OUTPUT}/
#mv ${SCRATCHDIR}/${SHORTNAME}.cons.fasta ${OUTPUT}

mv ${SCRATCHDIR}/* ${OUTPUT} && rm -rf ${SCRATCHDIR}



echo '******************************************'
echo "end at `date`"
echo '******************************************'
