#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1 # Un noeud par tache
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1  # Nb of threads we want to run on (il y a 32 CPU/noeud)
#SBATCH --partition=fast
#SBATCH --qos=fast
#SBATCH --export=ALL

OUTPUT='/home/palasser/projects/soutien_bioinfo/BAZILE/sequences_ABD_gene_RECQ4/genomes_diplo_tetra_hexa'
mkdir -p $OUTPUT
cd $OUTPUT
## DEMANDE: Recup des sequences des genes RECQ4 de REFSEQV1 suivant chez les genomes diplo, tetra, hexa: TraesCS2A02G304900 TraesCS2B02G321700 TraesCS2D02G303500

###########################################################################################################################################
# recupe des alignement gmap faits par Nathan
###########################################################################################################################################
for g in "Aet" "Tdur" "TraesTib" "Tru" "Zavv2" "ArinaLrFor" "Jagger" "Julius" "Lancer" "Landmark" "Mace" "Norin61" "Renan" "Spelt" "Stanley" "SY_Mattis";
do
    if [ $genome = "Aet" ];
        then
        gmap_bed=/home/napapon/results/gmapl_Aet_v4.0/IWGSC_v1.1_HC_20170706_cds.MergeAllChromD.9090.fasta-gmap-parsedOnlyBest_vs_Aet.bed

    elif [ $genome = "Tru" ];
        then
        gmap_bed=/home/napapon/results/gmapl_Trura_lastv/IWGSC_v1.1_HC_20170706_cds.MergeAllChromA.9090.fasta-gmap-parsedOnlyBest_vs_Tru_withTransloc7B4A.bed

    elif [ $genome = "Tdur" ];
        then
        gmap_bed=/home/napapon/results/gmapl_Tdurum/IWGSC_v1.1_HC_20170706_cds.MergeAllChromAB.9090.fasta-gmap-parsedOnlyBest_vs_Tdurum.bed

    elif [ $genome = "Zavv2" ];
        then
        gmap_bed=/home/napapon/results/gmapl_Zavv2.0/IWGSC_v1.1_HC_20170706_cds.MergeAllChrom.OneVariant.fasta-gmap-parsedOnlyBest_vs_Zavv2.gff3.bed

    elif [ $genome = "TraesTib" ];
        then    
        gmap_bed=/home/napapon/results/gmapl_TraesTib/IWGSC_v1.1_HC_20170706_cds.MergeAllChrom.OneVariant.fasta-gmap9090Cross-parsed_vs_TraesTib.gff3.bed
    
    elif [[ $genome = "ArinaLrFor" || $genome = "SY_Mattis" ]];
        then
        gmap_bed=/home/napapon/results/gmapl_${genome}/IWGSC_v1.1_HC_20170706_cds.MergeAllChrom_withTransloc7B5B.OneVariant.fasta-gmap9090Cross-parsed_vs_${genome}.gff3.bed

    else
        gmap_bed=/home/napapon/results/gmapl_${genome}/IWGSC_v1.1_HC_20170706_cds.MergeAllChrom.OneVariant.fasta-gmap9090Cross-parsed_vs_${genome}.gff3.bed
    fi

    egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' ${gmap_bed} |sed 's/ /\t/g' |cut -f1-4 |paste - <(grep $g $OUTPUT/path_genomes_diplo_tetra_hexa.tab |cut -f2) \
    >> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra_hexa.bed
done


egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Trura_lastv/IWGSC_v1.1_HC_20170706_cds.2A.OneVariant.fasta-gmap-parsedOnlyBest_vs_Trura_OnlyA2A_chr2A.gff3 \
> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_Triticum_urartu.gff3

##genome D Aegilops taushii
egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Aet_v4.0/IWGSC_v1.1_HC_20170706_cds.MergeAllChromD.9595.fasta-gmap-parsedOnlyBest_vs_Aet.bed \
|sed 's/ /\t/g' |cut -f1-4 |paste - <(echo "/home/napapon/Seq/Aegilops_tauschii/Aegilops_tauschii.Aet_v4.0.dna.allchromosome.fa") \
>> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed

egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Aet_v4.0/IWGSC_v1.1_HC_20170706_cds.2D.OneVariant.fasta-gmap-parsedOnlyBest_vs_Aet_OnlyD2D_2D.gff3 \
> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_Aegilops_tauschii.gff3

##genome AB Triticum dicoccoides
egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Zavv2.0/IWGSC_v1.1_HC_20170706_cds.MergeAllChrom.OneVariant.fasta-gmap-parsedOnlyBest_vs_Zavv2.gff3.bed \
|sed 's/ /\t/g' |cut -f1-4 |gawk '{print $0"\t/home/napapon/Seq/Triticum_dicoccoides/TrdicZavitanv2.0.fasta"}' \
>> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed

egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Zavv2.0/IWGSC_v1.1_HC_20170706_cds.2A.OneVariant.fasta-gmap-parsedOnlyBest_vs_Zavv2_2A_2A.gff3 \
> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_Triticum_dicoccoides.gff3

egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Zavv2.0/IWGSC_v1.1_HC_20170706_cds.2B.OneVariant.fasta-gmap-parsedOnlyBest_vs_Zavv2_2B_2B.gff3 \
>> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_Triticum_dicoccoides.gff3

##genome AB Triticum durum
egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Tdurum/IWGSC_v1.1_HC_20170706_cds.MergeAllChromAB.9090.fasta-gmap-parsedOnlyBest_vs_Tdurum.bed \
|sed 's/ /\t/g' |cut -f1-4 |gawk '{print $0"\t/home/napapon/Seq/Triticum_durum/Svevo.v1.0.april_2019.StdChrom.fna"}' \
>> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed

egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Tdurum/IWGSC_v1.1_HC_20170706_cds.2A.OneVariant.fasta-gmap-parsedOnlyBest_vs_TdurumChrom2A.gff3 \
> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_Triticum_durum.gff3

egrep 'TraesCS2A02G304900|TraesCS2B02G321700|TraesCS2D02G303500' /home/napapon/results/gmapl_Tdurum/IWGSC_v1.1_HC_20170706_cds.2B.OneVariant.fasta-gmap-parsedOnlyB2Best_vs_TdurumChrom2B.gff3 \
>> $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_Triticum_durum.gff3

###########################################################################################################################################
# recupe des sequences
###########################################################################################################################################
gawk '{print $0"\t"$5 }' $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed |cut -d'/' -f1-6,10 |sed -E 's/a\t\//a\t/' |sed 's/ID=//' |sed s'/lcl|//' \
> $OUTPUT/tmp && mv $OUTPUT/tmp $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed

module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 bioperl/1.7.0_rc5 bedtools/2.27.1
while read line;
do
    genome=$(echo $line |cut -d' ' -f6)
    path=$(echo $line |cut -d' ' -f5)
    gene_REFSEQV1=$(echo $line |cut -d' ' -f4)

    revcom.pl <(bedtools getfasta -name+ -bed <(grep $genome $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed |fgrep ${gene_REFSEQV1}) -fi $path) \
    > $OUTPUT/RECQ4_${genome}_ortho_de_${gene_REFSEQV1}_gene_startTOstop.fasta

done < $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed

ml gcc/8.1.0 cufflinks/2.2.1
sed -i s'/lcl|//' gmap_RECQ4_CDS.1_REFSEQV1_vs_Triticum_dicoccoides.gff3
while read line;
do
    genome=$(echo $line |cut -d' ' -f6)
    path=$(echo $line |cut -d' ' -f5)
    gene_REFSEQV1=$(echo $line |cut -d' ' -f4)

    gffread <(grep ${gene_REFSEQV1} $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_${genome}.gff3) -g ${path} -x $OUTPUT/RECQ4_${genome}_ortho_de_${gene_REFSEQV1}_cds.fasta -y $OUTPUT/RECQ4_${genome}_ortho_de_${gene_REFSEQV1}_pep.fasta

done < $OUTPUT/gmap_RECQ4_CDS.1_REFSEQV1_vs_genomes_diplo_tetra.bed
