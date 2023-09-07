#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1 # Un noeud par tache
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --partition=fast
#SBATCH --qos=fast
#SBATCH --export=ALL
#SBATCH --array=0-1

ml gcc/8.1.0 gmap_gsnap/v18.05.11 gdecTools/1.1 bioperl/1.7.0_rc5 cufflinks/2.2.1

OUTPUT='/home/palasser/projects/soutien_bioinfo/BAZILE/sequences_ABD_gene_RECQ4/genomes_ae_speltoides'
mkdir -p $OUTPUT
cd $OUTPUT
## DEMANDE: Recup des sequences des genes RECQ4 chez les genomes aegilops speltoides, aegilops mutica et triticum thimophevii dispo: REFSEQV1 TraesCS2A02G304900 TraesCS2B02G321700 TraesCS2D02G303500
## REFSEQV2: TraesCS2A03G0759500 TraesCS2B03G0837200 TraesCS2D03G0700700
##===> seul genome ae speltoides dispo pour deux accessions

genomes=("y2032" "GCA_021437245.1")
genome=${genomes[$SLURM_ARRAY_TASK_ID]}
fasta=("/storage/groups/gdec/shared/aegilops_speltoides/y2032/Y2032.updata.genome.fa" "/storage/groups/gdec/shared/aegilops_speltoides/GCA_021437245.1_speltoides/GCA_021437245.1_ASM2143724v1_genomic.fna")


gmapl -t 8 --cross-species -f 2 \
-D /storage/groups/gdec/shared/aegilops_speltoides/${genome}*/gmap_index \
-d gmap_index \
/home/palasser/projects/soutien_bioinfo/BAZILE/sequences_ABD_gene_RECQ4/RECQ4_chr2_IWGSC_refseqv2.1_annotation_200916_HC_cds.fasta \
> gmap_RECQ4_REFSEQV2_vs_${genome}.gff3

parseGmap.pl -c 70 -I 75 -r best -gmap gmap_RECQ4_REFSEQV2_vs_${genome}.gff3 \
-output gmap_RECQ4_REFSEQV2_vs_${genome}_parsed.gff3

gffread gmap_RECQ4_REFSEQV2_vs_${genome}_parsed.gff3 -g ${fasta[$SLURM_ARRAY_TASK_ID]} -x CDS_RECQ4_ortho_REFSEQV2_from_${genome}.fasta -y PEP_RECQ4_ortho_REFSEQV2_from_${genome}.fasta


### nouveaux genomes DOI: 10.1016/j.molp.2021.12.019 
#### https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA700474
## disponibles ici: /storage/groups/gdec/shared/aegilops_sitopsis