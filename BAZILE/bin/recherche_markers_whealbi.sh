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


