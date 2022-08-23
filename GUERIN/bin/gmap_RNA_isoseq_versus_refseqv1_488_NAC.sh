#!/bin/bash

#SBATCH --job-name=NACisoSEQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --partition=fast
#SBATCH --export=ALL
#SBATCH --array=0-5

module load gmap_gsnap/v18.05.11 triannotTools/1.2 gdecTools/1.1 bioperl/1.7.0_rc5


######## mapping evidence isoSeq   
#la premiere parenthese permet d'interpreter la sortie ls comme un array
ISO=($(\ls -1 /storage/groups/gdec/shared/triticum_aestivum/rnaseq/isoSeq-earlham2016/*.fna.gz |sed 's/.fna.gz//'))

sample=$(basename ${ISO[$SLURM_ARRAY_TASK_ID]})

gmap -t 8 --intronlength 5000 -f 2 -D /home/palasser/soutien_bioinfo/GUERIN/NAC_genes -d 488_NAC_CS_RefSeqv1_genes.gmapdb \
<(zcat /storage/groups/gdec/shared/triticum_aestivum/rnaseq/isoSeq-earlham2016/${sample}.fna.gz) > /home/palasser/soutien_bioinfo/GUERIN/NAC_genes/${sample}_isoseq_vs_refseqv1_488NAC.gff3

parseGmap.pl -c 25 -I 98 -r all -gmap /home/palasser/soutien_bioinfo/GUERIN/NAC_genes/${sample}_isoseq_vs_refseqv1_488NAC.gff3 \
-output /home/palasser/soutien_bioinfo/GUERIN/NAC_genes/${sample}_isoseq_vs_refseqv1_488NAC_parsed.gff3

