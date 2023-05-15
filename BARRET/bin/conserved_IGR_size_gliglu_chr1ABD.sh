#!/bin/bash
##SBATCH --time=1:00:00 #1h
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --array=0-2 
#SBATCH --cpus-per-task=8
#SBATCH --job-name=IGR_Ren
#SBATCH --partition=fast
#SBATCH --qos=fast

module load gmap_gsnap/v18.05.11 seqtk/1.3 ncbi-blast/2.11.0+

####################################################################################
# DEFINITION VARIABLES
####################################################################################
DATADIR=/storage/groups/gdec/shared/triticum_aestivum/wheatomics/renan/annot/v2/genes/v2021-09-08
OUTPUT=/home/palasser/projects/soutien_bioinfo/BARRET
cd $OUTPUT

##### regions d'interet
regions=("1A" "1B" "1D")
region=${regions[$SLURM_ARRAY_TASK_ID]}

# ####################################################################################
# EN AMONT
# ####################################################################################
# mkdir blastdb
# makeblastdb -parse_seqids -in cds_gliglu_${region}.fasta -dbtype nucl -out blastdb/cds_gliglu_${region}.blastdb

####################################################################################
# ANALYSES
####################################################################################

#gawk 'BEGIN { RS=">" } { if ($0~"TraesRN1") print RS $0}' $DATADIR/TaeRenan_refseqv2.1_genesHC_cds.fasta |sed -E 's/loc:chr.*[0-9]//' > TaeRenan_refseqv2.1_genesHC_cds.fasta

####################################################################################
bedr=/home/palasser/projects/storage_proteins_chr1/results/TaeRenan_refseq_v2.0/TaeRenan_refseq_v2.0_${region}S_region.bed
stop=$(cut -f3 $bedr)

egrep ^chr$region ${DATADIR}/TaeRenan_refseqv2.1_genesHC.gff3 |egrep $'\t''mRNA' |gawk -v s=$stop '{ if ($5<s) print $9}' |cut -d';' -f1 |sed 's/ID=//' > geneID_gliglu_${region}.txt

seqtk subseq TaeRenan_refseqv2.1_genesHC_cds.fasta geneID_gliglu_${region}.txt > cds_gliglu_${region}.fasta

## VERIF
wc -l geneID_gliglu_1*   =>  466 geneID_gliglu_1A.txt, 456 geneID_gliglu_1B.txt, 422 geneID_gliglu_1D.txt
grep -c '>' cds*fasta    =>  cds_gliglu_1A.fasta:466, cds_gliglu_1B.fasta:456, cds_gliglu_1D.fasta:422

for r in ${regions[@]/$region}; 
do 
    blastn -num_threads $SLURM_CPUS_PER_TASK -evalue 1e-3 -max_hsps 1 -max_target_seqs 1 \
    -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send evalue score qcovhsp" \
    -query cds_gliglu_${region}.fasta \
    -db blastdb/cds_gliglu_${r}.blastdb \
    -out Renan_cds_gliglu_${region}_VS_${r}.blastn

    gawk '{ if ($3>=90 && $14>=85) print $0 }' Renan_cds_gliglu_${region}_VS_${r}.blastn > tmp${region}_${r} && mv tmp${region}_${r} Renan_cds_gliglu_${region}_VS_${r}_parsed.blastn
done


# ####################################################################################
# ####################################################################################
join -t$'\t' -1 2 -2 1 <(cut -f1,2 Renan_cds_gliglu_1A_VS_1B_parsed2.blastn |sort -k2,2) <(cut -f1,2 Renan_cds_gliglu_1B_VS_1A_parsed2.blastn |sort -k1,1) |gawk '{if ($2==$3) print $1,$2 }' > Renan_cds_gliglu_1AB_1BA.txt

join -t$'\t' -1 2 -2 1 <(cut -f1,2 Renan_cds_gliglu_1A_VS_1D_parsed2.blastn |sort -k2,2) <(cut -f1,2 Renan_cds_gliglu_1D_VS_1A_parsed2.blastn |sort -k1,1) |gawk '{if ($2==$3) print $1,$2 }' > Renan_cds_gliglu_1AD_1DA.txt

join -t$'\t' -1 2 -2 2 <(sort -k2,2 Renan_cds_gliglu_1AB_1BA.txt) <(sort -k2,2 Renan_cds_gliglu_1AD_1DA.txt) > Renan_cds_gliglu_OG_1ABD.txt

sed -i 's/\.1//g' Renan_cds_gliglu_OG_1ABD.txt

####################################################################################
####################################################################################
join -t$'\t' -1 1 -2 4 <(sort -k1,1 Renan_cds_gliglu_OG_1ABD.txt) <(grep ^chr1A $DATADIR/TaeRenan_refseqv2.1_genesHC.bed |sort -k4,4) |gawk '{ print $1,$4,$5,$6,$2,$3}' > Renan_cds_gliglu_OG_1ABD_coordA.txt
join -t$'\t' -1 5 -2 4 <(sort -k5,5 Renan_cds_gliglu_OG_1ABD_coordA.txt) <(grep ^chr1B $DATADIR/TaeRenan_refseqv2.1_genesHC.bed |sort -k4,4) |gawk '{ print $2,$3,$4,$5,$1,$7,$8,$9,$6}' > Renan_cds_gliglu_OG_1ABD_coordAB.txt
join -t$'\t' -1 9 -2 4 <(sort -k9,9 Renan_cds_gliglu_OG_1ABD_coordAB.txt) <(grep ^chr1D $DATADIR/TaeRenan_refseqv2.1_genesHC.bed |sort -k4,4) |gawk '{ print $2,$3,$4,$5,$6,$7,$8,$9,$1,$10,$11,$12}' > Renan_cds_gliglu_OG_1ABD_coordABD.txt


sort -k3,3n Renan_cds_gliglu_OG_1ABD_coordABD.txt |gawk 'NR!=1 {e_before=a[4]} {split($0,a,FS)} { if (NR!=1) print $1,$3-e_before,$5,$6,$7,$8,$9,$10,$11,$12 ; if (NR==1) print $1,"-",$5,$6,$7,$8,$9,$10,$11,$12 }' > Renan_cds_gliglu_OG_1ABD_IGRsizeA.txt

sort -k5,5n Renan_cds_gliglu_OG_1ABD_IGRsizeA.txt |gawk 'NR!=1 {e_before=a[6]} {split($0,a,FS)} { if (NR!=1) print $1,$2,$3,$5-e_before,$7,$8,$9,$10 ; if (NR==1) print $1,$2,$3,"-",$7,$8,$9,$10 }' > Renan_cds_gliglu_OG_1ABD_IGRsizeAB.txt

sort -k7,7n Renan_cds_gliglu_OG_1ABD_IGRsizeAB.txt |gawk 'NR!=1 {e_before=a[8]} {split($0,a,FS)} { if (NR!=1) print $1,$2,$3,$4,$5,$7-e_before ; if (NR==1) print $1,$2,$3,$4,$5,"-" }' > Renan_cds_gliglu_OG_1ABD_IGRsizeABD.txtls
