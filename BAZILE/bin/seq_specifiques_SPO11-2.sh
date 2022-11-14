#!/bin/bash

### BUT de la manip, mail Jeanne BAZILE (13/09/2022): 
## Eteindre le gene SPO11-2 par transfo VIGS, sans eteindre SPO11-1.
# Le gene spo11-2 est encode TraesCS7B02G201200.1 and TraesCS7D02G296000.2, la copie A TraesCS7A02G300300.1, est non fonctionnelle
# cf Benyahya et al 2020 pour + d info: SPO11.2 is essential for programmed double-strand break formation during meiosis in bread wheat (Triticum aestivum L.).
# cahier des charges:
# 1 probe pour eteindre les 3 homeologues a la fois 
# probe de 100 a 200 bp max (taille max insert dasn virus), tendre vers 100bp
# probe ne doit pas contenir le site de l enzyme BSSHII
# probe doit etre specifique sur toute sa longueur (car dans le process decoupe en bout de 20bp environ), d'ou blastn word_size 5
# probe doit etre localisee dans un exon
# fournir 3 probes potentielles


mkdir /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11

fgrep 'Spo11' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.0/FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB

# SPO11-1: TraesCS5A02G391400 TraesCS5B02G396300 TraesCS5D02G401100
# SPO11-2: TraesCS7A02G300300 TraesCS7B02G201200 TraesCS7D02G296000
# SPOautre: TraesCS4A01G105600 TraesCS4B01G198400 TraesCS4D01G199200 TraesCS3D01G170600

#### SPO11 CDS
ml seqtk
rm /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11_refseq_v1.1_cds.fasta
for s in TraesCS5A02G391400 TraesCS5B02G396300 TraesCS5D02G401100 TraesCS7A02G300300 TraesCS7B02G201200 TraesCS7D02G296000 TraesCS4A02G105600 TraesCS4B02G198400 TraesCS4D02G199200 TraesCS3D02G170600;
do    
    gawk -v geneID=$s 'BEGIN { RS=">" } { if ($0~geneID) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_v1.1_HC_20170706_cds.fasta \
    |seqtk seq -l0 >> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11_refseq_v1.1_cds.fasta
done

ml ClustalOmega
clustalo -i /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11_refseq_v1.1_cds.fasta -t DNA --threads=8  --outfmt=clu \
> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11_refseq_v1.1_cds.clu


echo -e ">probe1_TraesCS7B02G201200.1\nGACGGGAAGCTTGTCACCCAGCGGGAGTTGTTCTACAAACTACTATCGGACTCACCCAAGTACTTCAGCTGTCAGCGCCATGTCAATCAAACCATCCAAG" \
> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2.fasta
echo -e ">probe2_TraesCS7B02G201200.1\nGAGGAGATGATCTGCAACTTATCTCTCAGAGTGCATTCCAAGAGTTGAAACCTCGTGATTTGCAGATCGCCAAAAGCTTGCTGTCATCTAAGTTTCTGCAGG" \
>> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2.fasta
echo -e ">probe3_TraesCS7B02G201200.1\nCCATTCTTGGACCTTCTGGGCATGCCATTACTGGGGATCTGAACCAATTGAGCAGATTAAATTTGTCTTCAGATGCCCGTTATCTCATTCTAGTGGAGAAGG" \
>> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2.fasta


ml ncbi-blast/2.11.0+
blastn -dust no -outfmt 6 -word_size 5 \
-query /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2.fasta \
-subject /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_v1.1_HC_20170706_mrna.fasta \
-out /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2_VS_IWGSC_v1.1_HC_mrna.blastn &

blastn -num_threads 8 -dust no -outfmt 6 -word_size 5 \
-query /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2.fasta \
-db /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta \
-out /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2_VS_IWGSC_v1.1.blastn

rm /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/gff_feature_blastn_match_probes_VIGS_vs_REFSEQV1.tab
while read line;
do
    probe=$(echo $line |cut -d' ' -f1)
    chrom=$(echo $line |cut -d' ' -f2)
    start=$(echo $line |cut -d' ' -f9)
    end=$(echo $line |cut -d' ' -f10)
    percentid=$(echo $line |cut -d' ' -f3)
    length=$(echo $line |cut -d' ' -f4)
    echo $probe" "$chrom"_"$start"_"$end" blastn results: "$percentid"%; "$length"bp" >> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/gff_feature_blastn_match_probes_VIGS_vs_REFSEQV1.tab
    gawk  -v chr=$chrom -v strt=$start -v nd=$end '{if ($1==chr && $4<strt && $5>nd) print}' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/IWGSC_v1.1_20170706.gff \
    |fgrep $'\t''gene' >> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/gff_feature_blastn_match_probes_VIGS_vs_REFSEQV1.tab
done < /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/probes_VIGS_spo11-2_VS_IWGSC_v1.1.blastn


#### SPO11-2 CDS
rm /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11-2_refseq_v1.1_cds.fasta
for s in TraesCS7A02G300300 TraesCS7B02G201200 TraesCS7D02G296000;
do    
    gawk -v geneID=$s 'BEGIN { RS=">" } { if ($0~geneID) print RS $0 }' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_v1.1_HC_20170706_cds.fasta \
    |seqtk seq -l0 >> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11-2_refseq_v1.1_cds.fasta
done

ml ClustalOmega
clustalo -i /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11-2_refseq_v1.1_cds.fasta -t DNA --threads=8  --outfmt=clu \
> /home/palasser/projects/soutien_bioinfo/BAZILE/SPO11/spo11-2_refseq_v1.1_cds.clu



### REFSEQV2
# egrep 'TraesCS5A02G391400\t|TraesCS5B02G396300\t|TraesCS5D02G401100\t|TraesCS7A02G300300\t|TraesCS7B02G201200\t|TraesCS7D02G296000\t|TraesCS4A02G105600\t|TraesCS4B02G198400\t|TraesCS4D02G199200\t|TraesCS3D02G170600'\
# /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV2/v2.1/annotation/genemodels_v200916/IWGSC_refseqv2.1_annotation_200916_IDmaping.csv

# egrep 'TraesCS7B03G0570200|TraesCS7D03G0707100|TraesCS7A03G0740200' ${DATABANK}/IWGSC_refseqv2.1_annotation_200916.gff3




