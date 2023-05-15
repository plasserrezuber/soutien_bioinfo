#!/bin/bash
#SBATCH --job-name=sbVIGS
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH -c 8
#SBATCH --array=0-17

ml samtools seqtk

##################################################################################
## WORKING DIRECTORIES
##################################################################################
DATADIR='/storage/groups/gdec/shared_mdc/WheatOmicsBam'
OUTPUT='/home/palasser/projects/soutien_bioinfo/MAGLIARASCHI'
mkdir -p ${OUTPUT}
cd ${OUTPUT}
mkdir bam

##################################################################################
## VARIABLES
## varietes d'interet = RECITAL et RENAN 
var='REN'
##################################################################################
SAMPLES=($(grep ^$var /home/palasser/projects/TEfusa/data/rnaseq_wheatomics/samples_list_MDC_RNAseq_WheatOMICS.txt |cut -f5 |tr '\n' ' '))
sample=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
file=($(\ls -1 ${DATADIR}/paired.postQC_CGC_${sample}*-Aligned.sortedByCoord.out.bam))
#########################################################

# fgrep 'TraesCS5A02G104900;' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/IWGSC_v1.1_20170706.gff \
# |grep $'\t''gene' |cut -f1,4,5,9 |cut -d';' -f1 |sed 's/ID=//' |gawk '{print $1,$2-1,$3,$4}' > VIGS_target_genes.bed
# fgrep 'TraesCS4A02G066200;' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/IWGSC_v1.1_20170706.gff \
# |grep $'\t''gene' |cut -f1,4,5,9 |cut -d';' -f1 |sed 's/ID=//' |gawk '{print $1,$2-1,$3,$4}' >> VIGS_target_genes.bed

# bedtools getfasta -nameOnly -fi /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta -bed VIGS_target_genes.bed -fo VIGS_target_genes.fasta

traitmt=$(grep $sample /home/palasser/projects/TEfusa/data/rnaseq_wheatomics/samples_list_MDC_RNAseq_WheatOMICS.txt |cut -f2)

samtools index -c $file

source /home/palasser/.bashrc
conda activate /home/palasser/apps/conda/envs/picard-2.27
##################################################################################
# extraire les reads mappes sur un chromosome pour les bam de plusieurs varietes
while read line;
do
    coord=$(gawk '{print $1":"$2"-"$3}' <(echo $line))
    gene=$(echo $line |cut -d' ' -f4)
    start=$(echo $line |cut -d' ' -f2)

    gawk -v strt=$start '{print $1,$2+strt,$3+strt,$4}' <(fgrep $gene VIGS_probes_coordinates_on_cDNA.bed) > ${gene}_VIGS_probes_coordinates_noIntron.bed

    grep $gene'[.;]' /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/annotation/v1.1/IWGSC_v1.1_20170706.gff \
    |grep $'\t''exon' |cut -f1,4,5  > ${gene}_cDNA.bed

    samtools view -@ ${SLURM_CPUS_PER_TASK} -bo bam/${gene}_${sample}_VIGS_target_${var}.bam $file $coord
    samtools index -c bam/${gene}_${sample}_VIGS_target_${var}.bam
    
    picard AddOrReplaceReadGroups \
    -I bam/${gene}_${sample}_VIGS_target_${var}.bam \
    -O bam/${gene}_${sample}_VIGS_target_${var}_GROUPS.bam \
    --RGID ${sample} \
    --RGLB ${var}_${traitmt} \
    --RGPL illumina \
    --RGPU ${var}_${traitmt} \
    --RGSM ${var}_${traitmt}

    samtools index -c bam/${gene}_${sample}_VIGS_target_${var}_GROUPS.bam

done < VIGS_target_genes.bed


# en dehors du job array
################################################################################################################
# source /home/palasser/.bashrc
# conda activate /home/palasser/apps/conda/envs/picard-2.27
# cd /home/palasser/projects/soutien_bioinfo/MAGLIARASCHI/bam
# for gene in $(cut -f4 /home/palasser/projects/soutien_bioinfo/MAGLIARASCHI/VIGS_target_genes.bed);
# do 
#     list=($(find ${gene}_*_VIGS_target_${var}_GROUPS.bam |sed 's/'$gene'/-I '$gene'/g'))
#     sbatch -p fast --qos=fast -c 2 --export=ALL --mem=64G --wrap="cd /home/palasser/projects/soutien_bioinfo/MAGLIARASCHI/bam; \
#     picard MergeSamFiles \
#     --USE_THREADING True \
#     $(echo ${list[@]}) \
#     -O merged_${gene}_VIGS_target_${var}_GROUPS.bam"
# done

# samtools index -c merged_TraesCS4A02G066200_VIGS_target_${var}_GROUPS.bam
# samtools index -c merged_TraesCS5A02G104900_VIGS_target_${var}_GROUPS.bam

# samtools faidx /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta chr4A > IWGSC_CSRefSeqv1_chr4A.fasta
# samtools faidx /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta chr5A > IWGSC_CSRefSeqv1_chr5A.fasta
# samtools faidx IWGSC_CSRefSeqv1_chr4A.fasta
# samtools faidx IWGSC_CSRefSeqv1_chr5A.fasta

ml gcc/8.1.0 bcftools/1.9
## si mapping avec bwa, mettre les options suivantes:
# -q 30: MAPQ reads >=30
# -Q 20: PHRED score >=20
# -C 50: adjust mapping quality; recommended:50, disable:0
# -m 10:  minimum number gapped reads of 10 for indel candidates

## -d: max per-file depth; avoids excessive memory usage [250] => sur RNAseq data, tres grande profondeur sur petites region, seuil pour ne pas saturer la RAM
bcftools mpileup  -R ../TraesCS4A02G066200_VIGS_probes_coordinates_withIntron.bed -o TraesCS4A02G066200_VIGS_target_${var}_mpileup.out -f /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta merged_TraesCS4A02G066200_VIGS_target_${var}_GROUPS.bam
bcftools mpileup  -R ../TraesCS5A02G104900_VIGS_probes_coordinates_withIntron.bed -o TraesCS5A02G104900_VIGS_target_${var}_mpileup.out -f /storage/groups/gdec/shared/triticum_aestivum/chinese_spring/iwgsc/REFSEQV1/fasta/IWGSC_CSRefSeqv1.fasta merged_TraesCS5A02G104900_VIGS_target_${var}_GROUPS.bam

bcftools call -cv TraesCS4A02G066200_VIGS_target_${var}_mpileup.out |grep -v 'ID=[1234MUH]' > TraesCS4A02G066200_VIGS_target_${var}.vcf
bcftools call -cv TraesCS5A02G104900_VIGS_target_${var}_mpileup.out |grep -v 'ID=[1234MUH]' > TraesCS5A02G104900_VIGS_target_${var}.vcf

ml bedtools
rm ../TraesCS4A02G066200_VIGS_target_${var}_vs_REFSEQV1.vcf
grep '#' TraesCS4A02G066200_VIGS_target_${var}.vcf > ../TraesCS4A02G066200_VIGS_target_${var}_vs_REFSEQV1.vcf
bedtools intersect -wa -wb -a ../TraesCS4A02G066200_VIGS_probes_coordinates_withIntron.bed -b TraesCS4A02G066200_VIGS_target_${var}.vcf \
|awk -v FS='\t' '{for (i=4; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' >> ../TraesCS4A02G066200_VIGS_target_${var}_vs_REFSEQV1.vcf

rm ../TraesCS5A02G104900_VIGS_target_${var}_vs_REFSEQV1.vcf
grep '#' TraesCS5A02G104900_VIGS_target_${var}.vcf > ../TraesCS5A02G104900_VIGS_target_${var}_vs_REFSEQV1.vcf
bedtools intersect -wa -wb -a ../TraesCS5A02G104900_VIGS_probes_coordinates_withIntron.bed -b TraesCS5A02G104900_VIGS_target_${var}.vcf \
|awk -v FS='\t' '{for (i=4; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' >> ../TraesCS5A02G104900_VIGS_target_${var}_vs_REFSEQV1.vcf