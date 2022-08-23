#!/bin/bash
#SBATCH --job-name=BLAST95
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --array=1-94

#############################################################
ml load ncbi-blast/2.11.0+ gdecTools/1.1  bioperl/1.7.0_rc5
#########################################################

OUTPUT='/home/newatil/blastn_pid95'
mkdir -p ${OUTPUT}
###################################################
##BLAST analysis
###################################################


for g in "FPM_i_1_A" "FPM_i_2_A" "FPM_m_1_D" "FPM_m_3_A" "FPM_m_3_D" "FPM_m_4_B" "FPM_m_4_D" "FPM_m_5_B" "FPM_m_5_D" "FPM_m_6_D" "FPM_m_7_D" "FPM_m_8_D" "HPM_1Ax" "HPM_1Ay" "HPM_1Bx"\
 "HPM_1By" "HPM_1Dx" "HPM_1Dy"
do
    FILE=$(\ls -1 /home/newatil/testPY/demux*_Glu_${g}_test.fasta |sed -n ${SLURM_ARRAY_TASK_ID}p)
    SHORTNAME=$(basename $FILE)

	blastn -num_threads 8 -dust no -max_hsps 1 -perc_identity 95.00 -strand plus -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
	-query /home/newatil/data/${g}_REFSEQV1.fasta \
	-subject ${FILE} \
	-out ${OUTPUT}/${SHORTNAME}_stplus.blastn

	blastn -num_threads 8 -dust no -max_hsps 1 -perc_identity 95.00 -strand minus -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
	-query /home/newatil/data/${g}_REFSEQV1.fasta \
	-subject ${FILE} \
	-out ${OUTPUT}/${SHORTNAME}_stminus.blastn
done
