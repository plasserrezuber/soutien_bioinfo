#!/bin/bash
#SBATCH --job-name=gampGlu
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=8

#####################
module load bioperl/1.7.0_rc5 gmap_gsnap/v18.05.11 bedtools samtools
##########################


OUTPUT='/home/palasser/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus'
mkdir ${OUTPUT}
cd ${OUTPUT}
#############################################

# for g in "hmw1by" "hmw1dy" "hmw1ax" "hmw1ay" "hmw1bx" "hmw1dx" "lmw1a_i_1" "lmw1a_i_2" "lmw1a_m_3" "lmw1b_m_4" "lmw1b_m_5" "lmw1d_m_1" "lmw1d_m_3" "lmw1d_m_4" "lmw1d_m_5" "lmw1d_m_6" "lmw1d_m_7" "lmw1d_m_8";
# do
#     sed -E 's/>Consensus\//\n>'$g' /' ${g}_cons.fasta >> 18genes_consensus_pbaa_revcom_clustlo_jalview.fasta
# done

#samtools faidx $OUTPUT/18genes_consensus_pbaa_revcom_clustlo_jalview.fasta
gmap_build -D $OUTPUT -d 18genes_consensus.gmapdb $OUTPUT/18genes_consensus_pbaa_revcom_clustlo_jalview.fasta

gmap --intronlength 500 -f 2 -D $OUTPUT -d 18genes_consensus.gmapdb \
/home/palasser/results/blast/gliglu_CHR1/REFSEQV2/IWGSC_refseq_v2.1_GLIGLU_CHR1_MANUALLY_CURATED_cds.fasta \
> $OUTPUT/gmap_18genes_CDS_gmap_vs_18consensus.gff3

parseGmap.pl -c 90 -I 85 -r best -gmap $OUTPUT/gmap_18genes_CDS_gmap_vs_18consensus.gff3 -output $OUTPUT/gmap_18genes_CDS_gmap_vs_18consensus_parse.gff3


