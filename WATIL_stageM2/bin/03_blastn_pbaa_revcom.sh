#!/bin/bash
#SBATCH --job-name=blst_pbaa
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --array=1-94

#############################################################
ml load ncbi-blast/2.11.0+ gdecTools/1.1  bioperl/1.7.0_rc5
#########################################################

OUTPUT='/home/newatil/blastn_pbaa'
mkdir -p ${OUTPUT}
cd ${OUTPUT}
###################################################
##BLAST analysis
#################################################

G=("HMW1Ax" "HMW1Ay" "HMW1Bx" "HMW1By" "HMW1Dx" "HMW1Dy" "LMW1A_i_1" "LMW1A_i_2" "LMW1A_m_3" "LMW1B_m_4" "LMW1B_m_5" "LMW1D_m_1" "LMW1D_m_3" "LMW1D_m_4" "LMW1D_m_5" "LMW1D_m_6" "LMW1D_m_7" "LMW1D_m_8")

i=0
for g in "hmw1ax" "hmw1bx" "hmw1dx" "hmw_1by" "hmw_1dy" "lmw1a_i_1" "lmw1a_i_2" "lmw1a_m_3" "lmw1b_m_4" "lmw1b_m_5" "lmw1d_m_1" "lmw1d_m_3" "lmw1d_m_4" "lmw1d_m_5" "lmw1d_m_6" "lmw1d_m_7" "lmw1d_m_8";
do 
  FILE=$(\ls -1 /home/newatil/multifa_split/sample*${G[i]}_cluster-0*fasta |sed -n ${SLURM_ARRAY_TASK_ID}p)
	SHORTNAME=$(basename ${FILE} |cut -d'.' -f1)
 
	blastn -num_threads 8 -dust no -max_hsps 1 -perc_identity 95.00 -strand plus \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
  -query /home/newatil/data/refseq/${g}_REFSEQV1.fasta \
  -subject ${FILE} \
  -out ${SHORTNAME}_stplus.blastn
  
  blastn -num_threads 8 -dust no -max_hsps 1 -perc_identity 95.00 -strand minus \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
  -query /home/newatil/data/refseq/${g}_REFSEQV1.fasta \
  -subject ${FILE} \
  -out ${SHORTNAME}_stminus.blastn
  
  ((i=i+1))
done

#cas particulier HMW1Ay
G="HMW1Ay"
cat /home/newatil/pbaa/pbaa_demux.*_passed_cluster_sequences.fasta |gawk 'BEGIN {RS=">" ; FS="\n"} {if ($1~"${g}") print ">"$0}' > /home/newatil/pbaa/${g}_assemblies_ALL.fasta

ml samtools
samtools faidx /home/newatil/pbaa/${g}_assemblies_ALL.fasta
gawk '{if ($2>6000) print $0"\t"$1}' /home/newatil/pbaa/${g}_assemblies_ALL.fasta.fai |cut -d'-' -f1-8 \
|awk -v OFS="\t" 'NR!=1 {before=a[6]} {split($0,a,FS)} { if ($6!=before) print}' |cut -f1 > /home/newatil/pbaa/${g}_assemblies_ALL_to_blast.txt

ml seqtk
seqtk subseq /home/newatil/pbaa/${g}_assemblies_ALL.fasta /home/newatil/pbaa/${g}_assemblies_ALL_to_blast.txt \
> /home/newatil/pbaa/${g}_assemblies_ALL_to_blast.fasta

rm /home/newatil/multifa_split/sample*HMW1Ay*.fasta
ml samtools
while read seq;
     do samtools faidx /home/newatil/pbaa/${g}_assemblies_ALL_to_blast.fasta $seq > /home/newatil/multifa_split/${seq}.fasta
 done < /home/newatil/pbaa/${g}_assemblies_ALL_to_blast.txt

FILE=$(\ls -1 /home/newatil/multifa_split/sample*${G}_*fasta |sed -n ${SLURM_ARRAY_TASK_ID}p)
SHORTNAME=$(basename ${FILE} |cut -d'.' -f1)

blastn -num_threads 8 -dust no -max_hsps 1 -perc_identity 95.00 -strand plus \
-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
-query /home/newatil/data/refseq/hmw1ay_REFSEQV1.fasta \
-subject ${FILE} \
-out ${SHORTNAME}_stplus.blastn

blastn -num_threads 8 -dust no -max_hsps 1 -perc_identity 95.00 -strand minus \
-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
-query /home/newatil/data/refseq/hmw1ay_REFSEQV1.fasta \
-subject ${FILE} \
-out ${SHORTNAME}_stminus.blastn


###################################################
## REVCOM for st minus vers /home/newatil/blastn_pbaa_rc 
###################################################
mkdir /home/newatil/blastn_pbaa_rc

module load gdecTools
files_minus=$(find /home/newatil/blastn_pbaa/sample*stminus.blastn -type f -not -empty |cut -d'.' -f1 |sed -n ${SLURM_ARRAY_TASK_ID}p)

for f in ${files_minus};
do
  seq=$(awk '{print $2}' ${f}.blastn)
  revcom.pl -o /home/newatil/blastn_pbaa_rc /home/newatil/multifa_split/${seq}.fasta
done

###################################################
## cp stplus vers /home/newatil/blastn_pbaa_rc 
###################################################

files_plus=$(find /home/newatil/blastn_pbaa/sample*stplus.blastn -type f -not -empty |cut -d'.' -f1 |sed 's/_stplus//' |sed -n ${SLURM_ARRAY_TASK_ID}p)
for f in ${files_plus};
do
  s=$(basename $f)
  cp /home/newatil/multifa_split/${s}.fasta /home/newatil/blastn_pbaa_rc
done


