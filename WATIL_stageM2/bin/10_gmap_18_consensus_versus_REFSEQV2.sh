#!/bin/bash
#SBATCH --job-name=Glugmap
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8
#SBATCH --array=0-17

#####################
module load gmap_gsnap/v18.05.11 bioperl/1.7.0_rc5 ncbi-blast/2.11.0+
##########################

GENOME='/storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1'
INPUT='/home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus'
OUTPUT='/home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18consensus_versus_REFSEQV2'
mkdir ${OUTPUT}
cd ${INPUT}


#############################################
GENES=("HMW1Ax" "HMW1Ay" "HMW1Bx" "HMW1By" "HMW1Dx" "HMW1Dy" "LMW1A_i_1" "LMW1A_i_2" "LMW1A_m_3" "LMW1B_m_4" "LMW1B_m_5" "LMW1D_m_1" "LMW1D_m_3" "LMW1D_m_4" "LMW1D_m_5" "LMW1D_m_6" "LMW1D_m_7" "LMW1D_m_8")
genes=("hmw1by" "hmw1dy" "hmw1ax" "hmw1ay" "hmw1bx" "hmw1dx" "lmw1a_i_1" "lmw1a_i_2" "lmw1a_m_3" "lmw1b_m_4" "lmw1b_m_5" "lmw1d_m_1" "lmw1d_m_3" "lmw1d_m_4" "lmw1d_m_5" "lmw1d_m_6" "lmw1d_m_7" "lmw1d_m_8")
g=${genes[$SLURM_ARRAY_TASK_ID]}

gmapl -t 8 --intronlength 300 -f 2 -D $GENOME -d gmap_index ${g}_cons.fasta > $OUTPUT/${g}_versus_REFSEQV2_gmap.gff3

parseGmap.pl -c 50 -I 90 -r best -gmap $OUTPUT/${g}_versus_REFSEQV2_gmap.gff3 -output $OUTPUT/${g}_versus_REFSEQV2_gmap_parsed.gff3


blastn -num_threads 8 -evalue 1e-10 -dust no -outfmt 6 \
-query ${g}_cons.fasta \
-db $GENOME/CS_pesudo_v2.1.fa \
-out $OUTPUT/${g}_versus_REFSEQV2.blastn

#######################################################################################
##### pour filtrer les fichiers de sortie Blastn:
#######################################################################################
# soit utiliser une fonction: function filterBlast () { awk -v id=$1 -v ov=$2 '{ if($3>=id && $4/$9>=ov){print} }'; } '$3>200 {print}
# "blastn xxxxxxxx | filterBlast 99 0.95 " pour filter les HSP sur 99percent id et 95percent overlap par exemple
gawk -F'\t' '{ if ($3>=90 && $4>400) {print $0} }' $OUTPUT/${g}_versus_REFSEQV2.blastn \
> $OUTPUT/${g}_versus_REFSEQV2_parsed.blastn