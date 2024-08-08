#!/bin/bash
#SBATCH --job-name=gampGlu
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1<

#####################
module load gcc/4.8.4 bioperl/1.7.0_rc5 gmap_gsnap/v18.05.11 samtools ncbi-blast/2.11.0+
##########################

########################################################
## ENVIRONNEMENT CONDA SNPEFF
########################################################
# conda create -p /home/palasser/apps/conda/envs/snpeff
# conda activate snpeff
# conda install -c bioconda snpeff=5.1 -y

########################################################
## verif si longueurs amplicons dans vcf identiques aux longueurs des consensus
########################################################

#grep '##contig' /home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18_x_94_msa2vcf/*.vcf
#cat /home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus/18genes_consensus.fasta.fai

##### pas pour tous les genes ==> produire le consensus emboss

INPUT=/home/palasser/projects/soutien_bioinfo/WATIL_stageM2
OUTPUT=/home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus
cd $OUTPUT

########################################################
## CONSENSUS SEQUENCES WITH EMBOSS
########################################################
ml java emboss/6.5.7

rm 18genes_consensus_emboss/18genes_emboss_cons.fasta
for g in "HMW1Ax" "HMW1Ay" "HMW1Bx" "HMW1By" "HMW1Dx" "HMW1Dy" "LMW1A_i_1" "LMW1A_i_2" "LMW1A_m_3" "LMW1B_m_4" "LMW1B_m_5" "LMW1D_m_1" "LMW1D_m_3" "LMW1D_m_4" "LMW1D_m_5" "LMW1D_m_6" "LMW1D_m_7" "LMW1D_m_8";
do
    showalign -show=d ${INPUT}/18_x_94_clustalo/pbaa_${g}_aling.fasta -outfile 18genes_consensus_emboss/${g}_aling.showalign
    cat <(echo ">"$g) <(grep Consensus 18genes_consensus_emboss/${g}_aling.showalign|sed 's/Consensus[[:space:]]*//') >> 18genes_consensus_emboss/18genes_emboss_cons.fasta
done

sed -iE 's/-/N/g' 18genes_consensus_emboss/18genes_emboss_cons.fasta
sed -iE 's/n/N/g' 18genes_consensus_emboss/18genes_emboss_cons.fasta

########################################################
# BLAST et GMAP CONSENSUS SEQUENCES VS _CDS_ RENAN  
########################################################

########################################################
blastn -dust no -max_hsps 1 -perc_identity 90.00 \
-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue score" \
-query /home/palasser/projects/storage_proteins_chr1/results/blast/gliglu_CHR1/TaeRenan_refseq_v2.0/TaeRenan_refseq_v2.0_GLIGLU_CHR1_MANUALLY_CURATED_cds.fasta \
-subject 18genes_consensus_emboss/18genes_emboss_cons.fasta \
-out 18genes_consensus_emboss/CDS_TaeRenan_VS_18genes_emboss_cons.blastn

########################################################
samtools faidx 18genes_consensus_emboss/18genes_emboss_cons.fasta
gmap_build -D 18genes_consensus_emboss -d 18genes_emboss_cons.gmapdb 18genes_consensus_emboss/18genes_emboss_cons.fasta

gmap --intronlength 500 -f 2 -D 18genes_consensus_emboss -d 18genes_emboss_cons.gmapdb \
/home/palasser/projects/storage_proteins_chr1/results/blast/gliglu_CHR1/TaeRenan_refseq_v2.0/TaeRenan_refseq_v2.0_GLIGLU_CHR1_MANUALLY_CURATED_cds.fasta \
> 18genes_consensus_emboss/gmap_CDS_TaeRenan_vs_18genes_emboss_cons.gff3

cd 18genes_consensus_emboss
parseGmap.pl -c 100 -I 90 -r best -gmap gmap_CDS_TaeRenan_vs_18genes_emboss_cons.gff3 -output gmap_CDS_TaeRenan_vs_18genes_emboss_cons_parsed.gff3

grep TraesRN1A01G00027600LC gmap_CDS_TaeRenan_vs_18genes_emboss_cons.gff3 |grep LMW1A_i_2 >> gmap_CDS_TaeRenan_vs_18genes_emboss_cons_parsed.gff3
grep TraesRN1D01G00036900 gmap_CDS_TaeRenan_vs_18genes_emboss_cons.gff3 |grep LMW1B_m_5 >> gmap_CDS_TaeRenan_vs_18genes_emboss_cons_parsed.gff3
grep -v TraesRN1A01G00034100 gmap_CDS_TaeRenan_vs_18genes_emboss_cons_parsed.gff3 > tmp && mv tmp gmap_CDS_TaeRenan_vs_18genes_emboss_cons_parsed.gff3


cd /home/palasser/projects/soutien_bioinfo/WATIL_stageM2

########################################################
# FASTA 18 SEQUENCES CONSENSUS for snpeff database building
########################################################
mkdir -p snpeff/data/18genes_consv0
cd snpeff

cp $OUTPUT/18genes_consensus_emboss/18genes_emboss_cons.fasta data/18genes_consv0/sequences.fa

########################################################
# GFF 18 SEQUENCES CONSENSUS for snpeff database building
########################################################
cut -d';' -f1 $OUTPUT/18genes_consensus_emboss/gmap_CDS_TaeRenan_vs_18genes_emboss_cons_parsed.gff3 > genes_tmp.gff
## reprise fichier pour n'avoir qu'une ligne pour exon et CDS
## boucle sur genes pour modif start stop des lignes gene et mRNA du gff (1 et longueur du gene grace au fichier 18genes_emboss_cons.fasta.fai)

rm data/18genes_consv0/genes.gff
while read line;
do
    g=$(echo $line |cut -d' ' -f1)
    l=$(echo $line |cut -d' ' -f2)

    gawk -v G=$g -v L=$l -v OFS='\t' '{ if ($1==G && $3=="gene") print G,$2,$3,$4,$5,$6,$7,$8,"gene_id "G"_GENE;"; if ($1==G && $3=="mRNA") print G,$2,$3,$4,$5,$6,$7,$8,"gene_id "G"_GENE; transcript_id "G"_CDS;"; if ($1==G && ($3=="exon" || $3=="CDS")) print G,$2,$3,$4,$5,$6,$7,".","transcript_id "G"_CDS;"}' genes_tmp.gff >> data/18genes_consv0/genes.gff
done < /home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus/18genes_consensus_emboss/18genes_emboss_cons.fasta.fai

sed -i 's/mRNA/transcript/' data/18genes_consv0/genes.gff
sed 's/d /d "/g' data/18genes_consv0/genes.gff |sed 's/;/";/g' > data/18genes_consv0/genes.gtf


#egrep -w "exon" data/18genes_consv0/genes.gff > tmp && mv tmp data/18genes_consv0/genes.gff

########################################################
# CHECK AND CREATE CDS FASTA FILE for snpeff database building
########################################################
 rm data/18genes_consv0/cds.fa
while read line;
do
    f=$(echo $line |cut -d' ' -f3)
    g=$(echo $line |cut -d' ' -f1)
    startcds=$(echo $line |cut -d' ' -f4)
    endcds=$(echo $line |cut -d' ' -f5)

    if [ $f = "CDS" ] 
    then
        samtools faidx ../18genes_CDS_gmap_vs_18consensus/18genes_consensus_emboss/18genes_emboss_cons.fasta $g:$startcds-$endcds |sed 's/:.*/_CDS/' >> data/18genes_consv0/cds.fa
    fi
done < data/18genes_consv0/genes.gff

#### RESULTAT OK !!:
###toutes les CDS obtenues a partir des consensus emboss ont bien un codon start en 5' si strand + dans data/18genes_consv0/genes.gff ou en 3' si strand - dans le gff3

########################################################
# CONFIG FILE for snpeff database building
########################################################
# touch 18genes_cons.config
    ## 18genes_cons genome, version 18genes_consv0
    ## 18genes_consv0.genome : 18genes_cons
    ## database.repository : /home/palasser/projects/soutien_bioinfo/WATIL_stageM2/snpeff

########################################################
# SNPEFF DATABASE DES 18 SEQUENCES CONSENSUS
########################################################

source ~/.bashrc
conda activate snpeff

#snpEff build -gff3 -c 18genes_cons.config 18genes_consv0 -noCheckProtein -v
snpEff build -gtf22 -c 18genes_cons.config 18genes_consv0 -noCheckProtein -noCheckCds -v

########################################################
# SNPEFF EFF
########################################################

for g in "HMW1Ax" "HMW1Ay" "HMW1Bx" "HMW1By" "HMW1Dx" "HMW1Dy" "LMW1A_i_1" "LMW1A_i_2" "LMW1A_m_3" "LMW1B_m_4" "LMW1B_m_5" "LMW1D_m_1" "LMW1D_m_3" "LMW1D_m_4" "LMW1D_m_5" "LMW1D_m_6" "LMW1D_m_7" "LMW1D_m_8";
do
    echo $g
    sed 's/chrUn/'$g'/' $INPUT/18_x_94_msa2vcf/pbaa_${g}_aling.fasta.vcf > ${g}_msa2vcf.vcf
    snpEff eff -c 18genes_cons.config -dataDir data 18genes_consv0 ${g}.vcf 1> ${g}_snpeff.vcf
    mv snpEff_genes.txt ${g}_snpEff_genes.txt
    mv snpEff_summary.html ${g}_snpEff_summary.html
    rm ${g}_msa2vcf.vcf
done

cat HMW1Ax_snpEff_genes.txt <(cat *_snpEff_genes.txt |grep -v '#' |grep -v 'HMW1Ax') > snpEff_genes.txt
 rm *_snpEff_genes.txt