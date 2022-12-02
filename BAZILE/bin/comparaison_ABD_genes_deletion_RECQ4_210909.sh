#!/bin/bash
module load bedtools/2.27.1

OUTPUT='/home/palasser/soutien_bioinfo/BAZILE'
mkdir $OUTPUT

## DEMANDE: Comparaison des genes chez Renan_pseudov2 entre les bornes de deletion suivantes
# Trouver le set commun des genes deletes entre les trois sous-genomes ABD
# Alignements des CDS des genes presents dans regions de deletion chez Renan de A contre B, A contre D et B contre D.
# Listes de genes absents pour 2, ou 3 regions de deletion


# fichier: $OUTPUT/ROI_Renan.txt
# chr2A   515429526       531057475
# chr2B   460671192       474547008
# chr2D   388551251       408812882

## regions sizes
# 15627949
# 13875816
# 20261631


cat /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/v2021-09-08/TaeRenan_refseqv2.0_genesHC_pep.fasta /storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/v2021-09-08/TaeRenan_refseqv2.0_genesLC_pep.fasta \
> /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseqv2.0_genesHC_LC_pep.fasta


## Obtention des seq proteiques
ml seqtk
seqtk subseq /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseqv2.0_genesHC_LC_pep.fasta \
<(cut -f1 liste_TaeRenan_refseq_v2.0_REFSEQV2_geneID_ROI_RECQ4.txt |grep '^TraesRN2A' |awk '{print $0".1"}') > TaeRenan_refseqv2.0_pep_ROI_2A_RECQ4.fasta
sed -i 's/\.1//' TaeRenan_refseqv2.0_pep_ROI_2A_RECQ4.fasta
sed -i 's/\./\*/g' TaeRenan_refseqv2.0_pep_ROI_2A_RECQ4.fasta

seqtk subseq /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseqv2.0_genesHC_LC_pep.fasta \
<(cut -f1 liste_TaeRenan_refseq_v2.0_REFSEQV2_geneID_ROI_RECQ4.txt |grep 'TraesRN2B' |awk '{print $0".1"}') > TaeRenan_refseqv2.0_pep_ROI_2B_RECQ4.fasta
sed -i 's/\.1//' TaeRenan_refseqv2.0_pep_ROI_2B_RECQ4.fasta
sed -i 's/\./\*/g' TaeRenan_refseqv2.0_pep_ROI_2B_RECQ4.fasta

seqtk subseq /home/palasser/data/RENAN_v2_pseudo/TaeRenan_refseqv2.0_genesHC_LC_pep.fasta \
<(cut -f1 liste_TaeRenan_refseq_v2.0_REFSEQV2_geneID_ROI_RECQ4.txt |grep 'TraesRN2D' |awk '{print $0".1"}') > TaeRenan_refseqv2.0_pep_ROI_2D_RECQ4.fasta
sed -i 's/\.1//' TaeRenan_refseqv2.0_pep_ROI_2D_RECQ4.fasta
sed -i 's/\./\*/g' TaeRenan_refseqv2.0_pep_ROI_2D_RECQ4.fasta


#### orthofinder
# source ~/.bashrc
# conda env list
# conda create -n orthofinder
# conda activate orthofinder
# conda search OrthoFinder
# conda config --show channels
# conda install orthofinder=2.5.4 -y  # option -y pour automatiquement valider les install
# orthofinder --help

conda activate orthofinder
sbatch -p fast -c 8 --export=ALL --wrap="ml ncbi-blast/2.11.0; orthofinder -a 8 -y -S blast -f /home/palasser/soutien_bioinfo/BAZILE/RECQ4_PEP_FASTA"

## mise en forme results orthofinder // BLASTP
sed -i 's/ //g' /home/palasser/soutien_bioinfo/BAZILE/RECQ4_PEP_FASTA/OrthoFinder/blastp_Nov02/Orthogroups/Orthogroups.tsv
gawk -v OFS='\t' '{print gsub(/Traes/,"Traes"),$0}' /home/palasser/soutien_bioinfo/BAZILE/RECQ4_PEP_FASTA/OrthoFinder/blastp_Nov02/Orthogroups/Orthogroups.tsv \
|gawk -v OFS='\t' '{print "_"gsub(/A/,"A",$0)gsub(/B/,"B",$0)gsub(/D/,"D",$0),$0}' |gawk '{print $3,$2,$1,$4,$5,$6}' |tr -d $'\r' |tr ' ' '\t' > $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv



##group 1 A ou B ou D
cat <(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,4 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,5 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,6 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,7 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,8 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,9 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,10 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,11 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,12 |gawk '{if ($4~"TraesRN2A") print $0}') \
<(cut -f1,2,3,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,13 |gawk '{if ($4~"TraesRN2A") print $0}') \
|sort -k2,2 > $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_A.tsv

##group 2 A ou B ou D
cat <(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,4 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,5 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,6 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,7 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,8 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,9 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,10 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,11 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,12 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,5 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,13 |gawk '{if ($4~"Traes") print $0}') |sort -k2,2 > $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_B.tsv

##group 3 A ou B ou D
cat <(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,4 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,5 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,6 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,7 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,8 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,9 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,10 |gawk '{if ($4~"Traes") print $0}') \
<(cut -f1,2,3,6 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS.tsv |tr ',' '\t' |cut -f1,2,3,11 |gawk '{if ($4~"Traes") print $0}') |sort -k2,2 > $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_D.tsv

cat $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_A.tsv $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_B.tsv $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_D.tsv |sort -k1,1 > $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_list.tsv

####RECUP NOMS GENES PUBLI POSSIBL SUR KNETMINER
sed -i 's/TRAES/Traes/' KnetMiner_Extract.txt
sed -i 's/02G/01G/' KnetMiner_Extract.txt

## join en deux fois pour gerer ordre des colonnes et comme colonnes d'annotation embetantes a parser
join -t$'\t' -a 1 -1 2 -2 1 <(join -t$'\t' -a 1 -1 2 -2 4 <(cut -f1,3 $OUTPUT/functional_annotation_REFSEQV1_TaeRenan_refseq_v2.0_ROI_RECQ4.txt |sort -k2,2) \
<(sort -k4,4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_list.tsv) |gawk -v OFS='\t' '{if ($3=="") {print $0,"NA","NA","NA"} else print $0}' |sed 's/\.1//' |sort -k2,2) \
<(cut -f1,2 $OUTPUT/KnetMiner_Extract.txt |sort -k1,1) |gawk -v OFS='\t' '{if ($6=="") {print $0,"NA"} else print $0}' > $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_KnetMiner_functionnal_annotation_tmp.tsv

join -t$'\t' -1 2 -2 1 <(sort -k2,2 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_KnetMiner_functionnal_annotation_tmp.tsv) \
<(cut -f3- $OUTPUT/functional_annotation_REFSEQV1_TaeRenan_refseq_v2.0_ROI_RECQ4.txt |sort -k1,1) \
> $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_KnetMiner_functionnal_annotation.tsv


# Stat rapides:
#cut -f1,2,4 TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_list.tsv |gawk '{ print $1,$2,substr($3,8,2) }' |sort |uniq -c |gawk '{print $1,$2,$3,$4}' > TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_list_ABD_grouped.tsv

#cut -f4 $OUTPUT/TaeRenan_refseqv2.0_RECQ4_HOMEOLOGS_KnetMiner_functionnal_annotation.tsv |sort -n |uniq -c
# 581 genes /834 seraient deletes pour 2 ou 3 sous-genomes a la fois
# 253 genes /834 n'appartiennent pas a un groupe de genes homeologues et ne seraient dons potentiellement pas deletes pour 2 ou 3 sous-genomes a la fois.

# Details
# nb genes         Nb genes/groupes
#     253                 0
#     172                 2
#     173                 3
#      44                  4
#      40                  5
#      51                  6
#       8                   8
#      18                  9
#      20                 10
#      22                 11
#      15                 15
#      18                 18
