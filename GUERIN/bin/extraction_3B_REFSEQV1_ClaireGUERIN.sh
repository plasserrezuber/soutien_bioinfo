#!/bin/bash
module load bedtools/2.27.1 samtools exonerate

# recuperation des sequences
bedtools getfasta -fi /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-bed /home/palasser/soutien_bioinfo/extraction_3B_REFSEQV1_ClaireGUERIN_210311.bed \
-fo /home/palasser/soutien_bioinfo/extraction_3B_REFSEQV1_ClaireGUERIN_210311.fasta

#verif taille attendue
samtools faidx /home/palasser/soutien_bioinfo/extraction_3B_REFSEQV1_ClaireGUERIN_210311.fasta

#explode 2 fata
fastaexplode /home/palasser/soutien_bioinfo/extraction_3B_REFSEQV1_ClaireGUERIN_210311.fasta -d /home/palasser/soutien_bioinfo