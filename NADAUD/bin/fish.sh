#!/bin/bash

# Faire une recherche de 2 sequences d'ADN g chez Renan (l'une fait 33 nucleotides et l'autre 50). Je t'explique : c'est pour mes manips d'hybridations
# in situ en fluorescence (GISH). La GISH permet de marquer entierement les chromosomes d'un meme genome (A en rouge, B en bleu, D en vert par exemple sur le ble tendre).
# En effet, je teste la faisabilite de la technique ND-FISH (FISH en conditions non denaturantes) en utilisant directement des oligos-sondes marques par un fluorochrome choisi et que l'on fait produire directement par invitrogen. Je me base sur un article de Tang et al 2018 qui donne un tableau de sequences d'oligo-sondes marques. Ils ont defini des sequences specifiques genome et specifiques chromosome et les ont teste sur Chinese spring. J'ai donc commande l'oligo B (qui permet de reperer tous les chrom du B) et l'oligo D (reconnait tous les chrom du D).
# Etant donne que les 1ers essais sont infructueux, je me pose la question de savoir si on retrouve bien cette sequence chez Renan (car j'utilise des lames de mitoses de Renan).
# Voici les sequences
# Oligo B :  5'  GGTTCAGGAATAGCCTCAGGAATTGGCTCAATT   3'
# Oligo D :  5'  TACGGGTGCCAAACGAGTGTCTGAAAGACTCCTCGAGAGGAAAATGCGAA  3'


##############################################################################################################################################
ml ncbi-blast/2.11.0+

blastn -num_threads 6 -dust no -perc_identity 70 -outfmt 6 \
-query /home/palasser/soutien_bioinfo/NADAUD/fish.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/blastdb/IWGSC_CSRefSeqv1 \
-out /home/palasser/soutien_bioinfo/NADAUD/fish_oligo_vs_Refseqv1.blastn &


blastn -num_threads 6 -dust no -perc_identity 70 -outfmt 6 \
-query /home/palasser/soutien_bioinfo/NADAUD/fish.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/CS_pesudo_v2.1.fa \
-out /home/palasser/soutien_bioinfo/NADAUD/fish_oligo_vs_Refseqv2.blastn &


blastn -num_threads 6 -dust no -perc_identity 70 -outfmt 6 \
-query /home/palasser/soutien_bioinfo/NADAUD/fish.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/blast/TaeRenan_refseq_v2.0.fa \
-out /home/palasser/soutien_bioinfo/NADAUD/fish_oligo_vs_RENAN_V2.blastn &

cut -f1-2 fish_oligo_vs_Refseqv1.blastn |sort |uniq -c
      1 Oligo_B chr3B
      1 Oligo_B chr5B
      1 Oligo_B chr6B
   1010 Oligo_D chr1D
   1544 Oligo_D chr2D
   1391 Oligo_D chr3D
   1102 Oligo_D chr4D
    988 Oligo_D chr5D
    786 Oligo_D chr6D
   1360 Oligo_D chr7D
    326 Oligo_D chrUn

cut -f1-2 fish_oligo_vs_Refseqv2.blastn |sort |uniq -c
      1 Oligo_B Chr3B
      1 Oligo_B Chr5B
      1 Oligo_B Chr6B
   1026 Oligo_D Chr1D
   1555 Oligo_D Chr2D
   1402 Oligo_D Chr3D
   1129 Oligo_D Chr4D
   1028 Oligo_D Chr5D
    890 Oligo_D Chr6D
   1392 Oligo_D Chr7D
    140 Oligo_D ChrUnknown


cut -f1-2 fish_oligo_vs_RENAN_V2.blastn |sort |uniq -c
      1 Oligo_B chr3B
      1 Oligo_B chr5B
      1 Oligo_B chr6B
      6 Oligo_D chr1B
    992 Oligo_D chr1D
   1429 Oligo_D chr2D
      1 Oligo_D chr3A
   1397 Oligo_D chr3D
   1139 Oligo_D chr4D
   1075 Oligo_D chr5D
    847 Oligo_D chr6D
   1287 Oligo_D chr7D
