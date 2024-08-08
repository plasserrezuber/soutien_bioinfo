#!/bin/bash
#SBATCH --job-name=gampGlu
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1

#####################
module load gcc/4.8.4 samtools/1.3
##########################

# sed 's/::.*//' 18genes_consensus_pbaa_revcom_clustlo_jalview.fasta > 18genes_consensus.fasta
# samtools faidx 18genes_consensus.fasta

OUTPUT='/home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus/18genes_consensus_jalview'
cd ${OUTPUT}
#############################################
## CONSENSUS JALVIEW
#############################################

#GENES=("HMW1Ax" "HMW1Ay" "HMW1Bx" "HMW1By" "HMW1Dx" "HMW1Dy" "LMW1A_i_1" "LMW1A_i_2" "LMW1A_m_3" "LMW1B_m_4" "LMW1B_m_5" "LMW1D_m_1" "LMW1D_m_3" "LMW1D_m_4" "LMW1D_m_5" "LMW1D_m_6" "LMW1D_m_7" "LMW1D_m_8")

while read line;
do
    g=$(echo $line |cut -d' ' -f1 |sed -E 's/_[[:digit:]]*$//' |tr '[:upper:]' '[:lower:]')

    snp=$(echo $line |cut -d' ' -f1)

    G=$(echo $line |cut -d' ' -f1 |sed -E 's/_[[:digit:]]*$//')
    pos=$(echo $line |cut -d' ' -f1 |sed 's/'$G'_//')

    REF=$(echo $line |cut -d' ' -f2)
    ALT=$(echo $line |cut -d' ' -f3)

    if [ $pos -gt 100 ] && [ ${#REF} -eq 1 ] && [ ${#ALT} -eq 1 ] ;
    then
        let "start5=$pos - 101"
        let "end5=$pos - 1"
        let "start3=$pos + 1"
        let "end3=$pos + 101"

        seq5=$(samtools faidx 18genes_consensus.fasta $g:$start5-$end5 |grep -v '>' |tr -d '\n')
        seq3=$(samtools faidx 18genes_consensus.fasta $g:$start3-$end3 |grep -v '>' |tr -d '\n')

        echo -e $snp"\t"$seq5"["$REF"/"$ALT"]"$seq3 >> glutenin_SNP_context_seq.tsv
    fi
done < <(grep -v 'marker' ../code_markers_18_genes_glu_FSOV_Exige.tab)


#############################################
## CONSENSUS EMBOSS
#############################################

OUTPUT='/home/palasser/projects/soutien_bioinfo/WATIL_stageM2/18genes_CDS_gmap_vs_18consensus/18genes_consensus_emboss'
cd ${OUTPUT}

while read line;
do
    G=$(echo $line |cut -d' ' -f1 |sed -E 's/_[[:digit:]]*$//')

    snp=$(echo $line |cut -d' ' -f1)

    pos=$(echo $line |cut -d' ' -f1 |sed 's/'$G'_//')

    REF=$(echo $line |cut -d' ' -f2)
    ALT=$(echo $line |cut -d' ' -f3)

    if [ $pos -gt 100 ] && [ ${#REF} -eq 1 ] && [ ${#ALT} -eq 1 ] ;
    then
        let "start5=$pos - 101"
        let "end5=$pos - 1"
        let "start3=$pos + 1"
        let "end3=$pos + 101"

        seq5=$(samtools faidx 18genes_emboss_cons.fasta $G:$start5-$end5 |grep -v '>' |tr -d '\n')
        seq3=$(samtools faidx 18genes_emboss_cons.fasta $G:$start3-$end3 |grep -v '>' |tr -d '\n')

        echo -e $snp"\t"$seq5"["$REF"/"$ALT"]"$seq3 >> glutenin_SNP_context_seq_EMBOSS.tsv
    fi
done < <(grep -v 'marker' ../code_markers_18_genes_glu_FSOV_Exige.tab)