#!/bin/bash

######################## activate conda environnement #########
conda activate msa2vcf


####################### input files ###########################
###fichiers input:
files=$(ls -1 /home/newatil/pbaa_rc_clustalo/pbaa*aling.fasta)


##########################################################
########## generate vcf files with msa2vcf ###############
##########################################################
## script execute en local car msa2vcf non installe sur hpc2
for f in ${files};
do 
  java -jar dist/msa2vcf.jar ${f} > ${f}.vcf
done


