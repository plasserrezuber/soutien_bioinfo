#!/bin/bash
OUTPUT='/home/palasser/soutien_bioinfo/GEORGES'
mkdir $OUTPUT

ml seqtk ncbi-blast/2.11.0+ 

###########################################################################################################################################
#### Join REFSEQV1 REFSEQV2: ## recup des geneID REFSEQV2 a partir de la liste de geneID trouves sur KnetMiner
###########################################################################################################################################
join -t$'\t' -1 1 -2 1 <(cut -f1 $OUTPUT/liste_geneID_Msh.txt |sort) \
<(sort -k1,1 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_IDmaping.csv) \
> $OUTPUT/correspondance_geneID_Msh_REFSEQV1_REFSEQV2.txt

seqtk subseq /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/genes_2017July06/IWGSC_v1.1_HC_20170706_transcripts.fasta \
<(cut -f1 $OUTPUT/liste_geneID_Msh.txt |gawk '{print $0".1"}') > sequences_Msh_REFSEQV1.fasta

seqtk subseq /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC_mrna.fasta \
<(cut -f2 $OUTPUT/correspondance_geneID_Msh_REFSEQV1_REFSEQV2.txt |gawk '{print $0".1"}') > sequences_Msh_REFSEQV2.fasta


###########################################################################################################################################
#### alignement blastn Arabido th versus REFSEQV2
###########################################################################################################################################
blastp -num_threads 8 -evalue 1e-5 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send evalue score" \
-query $OUTPUT/seq_pep_Msh_Ath.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_HC_pep.valid.fasta \
-out $OUTPUT/pep_Msh_Ath_vs_REFSEQV2.blastp

echo -e "query_REFSEQV2\tquery_REFSEQV1\tsubject\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tevalue\tscore" > pep_Msh_Ath_vs_REFSEQV2_corresp_REFSEQV1.blastp
join -t$'\t' -1 2 -2 2 <(sort -k 2,2 /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_IDmaping.csv |gawk '{print $0".1"}') \
<(sort -k 2,2 $OUTPUT/pep_Msh_Ath_vs_REFSEQV2.blastp) |sort -k3,3 -k14,14nr |gawk -v OFS='\t' '{if ($4>=50) {print $0}}' >> pep_Msh_Ath_vs_REFSEQV2_corresp_REFSEQV1.blastp