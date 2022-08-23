#!/bin/bash
##SBATCH --time=1:00:00 #1h
#SBATCH --nodes=1 # Un noeud par tache
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --array=0-9 # chaque mapping gmap va etre execute sur 1 noeud (10 noeuds en parallel)
#SBATCH --cpus-per-task=8  # Nb of threads we want to run on (il y a 32 CPU/noeud)
#SBATCH --job-name=bedtools_getfasta
#SBATCH --partition=fast

module restore

blastdbcmd -db $TAE/chinese_spring/iwgsc/REFSEQV1/blastdb/IWGSC_v1.1_mrna -entry_batch \
<(cut -f1 $HOME/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.TAB | awk '{print $0".1"}') > $HOME/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.fna

####################################################################################
##### genomes d'interet
genomes=("arinaLrFor" "jagger" "julius" "Lancer" "landmark" "mace" "Norin61" "spet" "stanley" "SY_mattis")
# genome au singulier va prendre la valeur contenue dans l'array "genomes" selon l'indice de la var d'envir $SLURM_ARRAY_TASK_ID
genome=${genomes[$SLURM_ARRAY_TASK_ID]}

####################################################################################
##### gmapl
module load gmap_gsnap/v18.05.11
# option --intronlength 5000 pour eviter les hit de 40kb....
# option -f 2 pour format de sortie type ggf3_gene

mkdir $HOME/results/gmap/${genome}
gmapl --intronlength 5000 -f 2 -D $TAE/10wheatgenomes/${genome}/${genome}.gmapdb \
-d ${genome}.gmapdb $HOME/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.fna > $HOME/results/gmap/${genome}/gmap_gli_glu_refseqv1.0_FuncAnnot_sur_${genome}.gff3

####################################################################################
##### parseGmap et mise en forme
# fichier de sortie gmap parse avec couverture >=90% et Identity >=75% et on garde le meilleur hit pour chaque query après ces deux filtres
parseGmap.pl -gmap $HOME/results/gmap/${genome}/gmap_gli_glu_refseqv1.0_FuncAnnot_sur_${genome}.gff3 \
-output $HOME/results/gmap/${genome}/gmap_gli_glu_refseqv1.0_FuncAnnot_sur_${genome}_parse.gff3 -c 90 -I 75 -r best

# pour ne garder que les "path" avec info du hit (lignes avec feature=mRNA)
gawk -F'\t' '{print $1"\t"$4-1"\t"$5"\t"$9}' <(grep 'mRNA' $HOME/results/gmap/${genome}/gmap_gli_glu_refseqv1.0_FuncAnnot_sur_${genome}_parse.gff3) \
| cut -d ';' -f1 | sed "s/ID\=//" | sed "s/.mrna1//" > $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}.bed

# recup de l'info d'annotation avec un join entre liste de travail et le bed obtenu
join -1 4 -2 1 <(sort -k4,4 $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}.bed | sed "s/lcl|//") \
<(sort -k1,1 $HOME/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.TAB | cut -f1,4 | awk '{print $1".1\t" $2}') \
| awk 'BEGIN{FS=" ";OFS="\t"} {print $2,$3,$4,$1" "$5}' > $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}_annot.bed

# pour trier par ordre croissant chrom et start avant utilisation bedtools
sort -V -k1,1 -k2,2 $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}_annot.bed \
> temp && mv temp $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}_annot.bed

####################################################################################
##### betools getfasta
module load bedtools/2.27.1
# recuperation des fasta des genes gli_glu pour les 10genomes a fournir a Mireille
mkdir $HOME/results/bedtools/tenGenomes

bedtools slop -l 6000 -r 6000 -i $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}_annot.bed -g <(cut -f1,2 $CS/fasta/IWGSC_CSRefSeqv1.fasta.fai) \
> $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}_annot_6k.bed

bedtools getfasta -name+ -fi $TAE/10wheatgenomes/${genome}/*.fasta -bed $HOME/results/gmap/${genome}/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}_annot_6k.bed \
-fo $HOME/results/bedtools/tenGenomes/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_${genome}.fna

# recup sequences Chinese-spring:
# bedtools getfasta -name+ -fi $CS/fasta/IWGSC_CSRefSeqv1.fasta -bed gli_glu_refseqv1.0_FuncAnnot_CS_annot_6k.bed -fo $HOME/results/bedtools/tenGenomes/gli_glu_refseqv1.0_FuncAnnot_gmap_sur_CS.fna

# wc -l $HOME/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.TAB
# 161
# grep -c '>' $HOME/data/gli_glu_iwgsc_refseqv1.0_FuncAnnot_v1.1_completed_v1.0.fna
# 143

# grep -c '>' *.fna
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_arinaLrFor.fna:113
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_jagger.fna:107
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_julius.fna:107
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_Lancer.fna:118
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_landmark.fna:116
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_mace.fna:114
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_Norin61.fna:111
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_spet.fna:114
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_stanley.fna:111
# gli_glu_refseqv1.0_FuncAnnot_gmap_sur_SY_mattis.fna:118
