#!/bin/bash

### BUT de la manip: extraire les variations obervees autour du gene RECQ4 chez un panel d'accession pour etablir les haplotypes et caracteriser la variabilite nucleotidique du gene RECQ4 et de sa region

## liste des genes avant/apres RECQ4 pur REFSEQV1_v1.1: recherche des Traes des genes voisins de RECQ4 car a part pour le chr2D, RECQ4 n'est pas retrouve dans le vcf de whealbi
egrep -C 60 'TraesCS2A02G304900;|TraesCS2B02G321700;|TraesCS2D02G303500;' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/IWGSC_v1.1_20170706.gff |fgrep -w gene

## genes les plus proches de RECQ4 (gene voisin direct ou RECQ4 directement pour chr2D) retenu pour l'extraction de 25 variations avant/apres du vcf whealbi
egrep 'CLOSEST_HC_GENE=TraesCS2A01G304800|CLOSEST_HC_GENE=TraesCS2B01G321600|CLOSEST_HC_GENE=TraesCS2D01G303500' \
/storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/markers/iwgsc_refseqv1.0_Whealbi_GWAS/whealbi.imputed.vcf |cut -f1-8

### extraction des 25 marqueurs avaat/apres ces genes, donc environnant RECQ4
grep '#' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/markers/iwgsc_refseqv1.0_Whealbi_GWAS/whealbi.imputed.vcf \
> /home/palasser/soutien_bioinfo/BAZILE/extraction_around_RECQ4_FROM_whealbi.imputed.vcf
egrep -C 25 'CLOSEST_HC_GENE=TraesCS2A01G304800|CLOSEST_HC_GENE=TraesCS2B01G321600|CLOSEST_HC_GENE=TraesCS2D01G303500' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/markers/iwgsc_refseqv1.0_Whealbi_GWAS/whealbi.imputed.vcf \
>> /home/palasser/soutien_bioinfo/BAZILE/extraction_around_RECQ4_FROM_whealbi.imputed.vcf




#### stat
egrep -C 25 'CLOSEST_HC_GENE=TraesCS2A01G304800|CLOSEST_HC_GENE=TraesCS2B01G321600|CLOSEST_HC_GENE=TraesCS2D01G303500' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/markers/iwgsc_refseqv1.0_Whealbi_GWAS/whealbi.imputed.vcf |fgrep -w chr2A |wc -l
64
egrep -C 25 'CLOSEST_HC_GENE=TraesCS2A01G304800|CLOSEST_HC_GENE=TraesCS2B01G321600|CLOSEST_HC_GENE=TraesCS2D01G303500' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/markers/iwgsc_refseqv1.0_Whealbi_GWAS/whealbi.imputed.vcf |fgrep -w chr2B |wc -l
59
egrep -C 25 'CLOSEST_HC_GENE=TraesCS2A01G304800|CLOSEST_HC_GENE=TraesCS2B01G321600|CLOSEST_HC_GENE=TraesCS2D01G303500' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/markers/iwgsc_refseqv1.0_Whealbi_GWAS/whealbi.imputed.vcf |fgrep -w chr2D |wc -l
55


### position variation 2B sur REFSEQV1
ml ncbi-blast/2.11.0+
blastn -num_threads 8 -dust no -outfmt 6 -word_size 5 \
-query /home/palasser/projects/soutien_bioinfo/BAZILE/sequences_ABD_gene_RECQ4/variation_RECQ4_2B.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-out /home/palasser/projects/soutien_bioinfo/BAZILE/sequences_ABD_gene_RECQ4/variation_RECQ4_2B.blastn

# position su snp dans geneTraesCS2B02G321700 chr2B:458849313..458858559 (- strand)  = chr2B:458854556
# Info snp:Project: http://wheatgenomics.plantpath.ksu.edu/1000EC/ 
# Name: chr2B_scaffold13731-2_13318540
# Type: SNV
# Description: SNV C -> T
# Position: chr2B:458854556..458854556
# Length: 1 bp
# Attributes: AC 32 AN 1620
# alternative_alleles: T
# description: SNV C -> T
# reference_allele: C
# seq_id: chr2B

#fichiers fourni a Jeanne Bazile:
/home/palasser/projects/soutien_bioinfo/BAZILE/extraction_snp_RECQ4_2B.vcf
http://wheatgenomics.plantpath.ksu.edu/1000EC/files/PassportData_160809.xlsx