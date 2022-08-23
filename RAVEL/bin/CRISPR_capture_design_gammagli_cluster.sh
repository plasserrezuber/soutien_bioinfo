#!/bin/bash
##SBATCH --time=1:00:00 #1h
#SBATCH --nodes=1 # Un noeud par tache
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8  # Nb of threads we want to run on (il y a 32 CPU/noeud)
#SBATCH --job-name=crispgam
#SBATCH --partition=fast

module load bedtools/2.27.1 ncbi-blast/2.11.0+

DATABANK1='/storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1'
DATABANK2='/home/palasser/data/RENAN_v2_pseudo/'

OUTPUT='/home/palasser/soutien_bioinfo/RAVEL'

###########################################################################################################################################
bedtools getfasta -name -fi $DATABANK1/CS_pesudo_v2.1.fa -bed $OUTPUT/CS_REFSEQV2_gammagli_seq_for_capture_design.bed \
-fo $OUTPUT/CS_REFSEQV2_gammagli_seq_for_capture_design.fasta

###########################################################################################################################################
bedtools getfasta -name -fi $DATABANK2/TaeRenan_refseq_v2.0.fa -bed $OUTPUT/TaeRenan_refseq_v2.0_gammagli_seq_for_capture_design.bed \
-fo $OUTPUT/TaeRenan_refseq_v2.0_gammagli_seq_for_capture_design.fasta


### BLAST CS
###########################################################################################################################################
blastn -num_threads 8 -evalue 1e-10 -perc_identity 90 -outfmt 0 \
-query $OUTPUT/CS_REFSEQV2_gammagli_seq_for_capture_design.fasta \
-db $DATABANK1/CS_pesudo_v2.1.fa \
-out $OUTPUT/CS_REFSEQV2_gammagli_capture_design_vs_REFSEQV2.blastn

###########################################################################################################################################
blastn -num_threads 8 -evalue 1e-10 -outfmt 6 \
-query $OUTPUT/CS_REFSEQV2_gammagli_seq_for_capture_design.fasta \
-db $DATABANK1/CS_pesudo_v2.1.fa \
-out $OUTPUT/CS_REFSEQV2_gammagli_capture_design_vs_REFSEQV2_tab.blastn

### BLAST RENAN
###########################################################################################################################################
blastn -num_threads 8 -evalue 1e-10 -perc_identity 90 -outfmt 0 \
-query $OUTPUT/TaeRenan_refseq_v2.0_gammagli_seq_for_capture_design.fasta \
-db $DATABANK2/blastdb/TaeRenan_refseq_v2.0.fa \
-out $OUTPUT/TaeRenan_refseq_v2.0_gammagli_capture_design_vs_TaeRenanv2.blastn

###########################################################################################################################################
blastn -num_threads 8 -evalue 1e-10 -outfmt 6 \
-query $OUTPUT/TaeRenan_refseq_v2.0_gammagli_seq_for_capture_design.fasta \
-db $DATABANK2/blastdb/TaeRenan_refseq_v2.0.fa \
-out $OUTPUT/TaeRenan_refseq_v2.0_gammagli_capture_design_vs_TaeRenanv2_tab.blastn