#!/bin/bash

#### jvarkit msa2vcf 
# source ~/.bashrc
# conda env list
# conda config --show channels
# conda create -p /home/palasser/apps/conda/envs/msa2vcf
# conda install -c bioconda jvarkit-msa2vcf

# cd /home/palasser/apps/conda/envs/msa2vcf/bin
# for f in /home/palasser/apps/conda/envs/msa2vcf/share/jvarkit-msa2vcf-201904251722-1/msa2vcf*; do ln -s $f; done
# java -jar msa2vcf.jar --help
# conda deactivate

######################## activate conda environnement #########
conda activate msa2vcf


####################### input files ###########################
###fichiers input:
files=$(ls -1 /home/newatil/pbaa_rc_clustalo/pbaa*aling.fasta)


##########################################################
########## generate vcf files with msa2vcf ###############
##########################################################
## script execute en local sur PC Nezha Watil.....
#NB: sans option de consensus, msa2vcf considere comme chrUn la longueur totale de l'alignement multiple
# les coordonnees sont celles sur le chrUn
# l'allele REF est l'allele majoritaire, les alleles ALT sont les alleles minoritaires

for f in ${files};
do 
  java -jar dist/msa2vcf.jar ${f} > ${f}.vcf
done

