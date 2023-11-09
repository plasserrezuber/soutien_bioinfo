#!/bin/bash

#SBATCH --job-name=archivage
#SBATCH -o archivage-%A_%a_task.out
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --partition=normal
#SBATCH --export=all

if [ "$#" -ne 1 ]; then
        echo -e "Need 1 argument :\n"
        echo -e "-> The path of the files to archive without the   /   at the end : eg. 84087_20230613_141110/1_A01 \n"
        exit
fi

echo -e "running on `hostname`\n"

#SetUpFolders


workingdir='/storage/scratch/'$SLURM_JOBID

mkdir $workingdir && cd $workingdir 

echo -e "The workingdir's path is :" $(pwd) "\n"


#LoadModules

module load python s3cmd goofys

mypath=${1?Usage: $0 <mypath>}

#Set up the mount



mkdir pacbio-data-provisoire

module load python goofys

export AWS_ACCESS_KEY_ID='b27b3d0a320941e3b4afab555faf0f2f'
export AWS_SECRET_ACCESS_KEY='f323f20b0dfd436295c44c0ded78dd99'

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire ./pacbio-data-provisoire


#Set up the name of the archive

archive=$(fgrep -m1 'ConsensusReadSet' pacbio-data-provisoire/$mypath/pb_formats/*.hifi_reads.consensusreadset.xml |cut -d' ' -f4 | cut -d'"' -f2| sed 's/-Cell[1-9]//g')
smrtcell_name=$(fgrep -m1 'ResourceId' pacbio-data-provisoire/$mypath/pb_formats/*.hifi_reads.consensusreadset.xml | cut -d' ' -f4 | cut -d/ -f3 | cut -d'.' -f1)


## utilisation de lien symboliques du point de montage ./pacbio-data-provisoire vers le workingdir (on a fait un cd en débit de script) pour éviter une écriture directe dans le point de montage
# le tar.gz est ecrit en dehors du point de montage, puis verse avec s3cmd put sur oscar
# (qui serait sauvegardée sur oscar également, ce qui soliciterai bcp d'appels reseau + des snapshots de sauvegarde)
for i in pacbio-data-provisoire/$mypath/*/ ; do ln -sr $i; done

tar --exclude='*.hifi_reads.bc*.bam' --exclude='*.hifi_reads.default.bam' -czvhf  ${smrtcell_name}.${archive}.metadata.tar.gz fail_reads/ hifi_reads/ metadata/ pb_formats/ statistics/

echo -e "\n"

s3cmd -c ~/.s3cfg_gentyane put ${smrtcell_name}.${archive}.metadata.tar.gz s3://pacbio-data-provisoire/${mypath}/${smrtcell_name}.${archive}.metadata.tar.gz

echo -e "\n"

echo -e "##############################################################################"
echo -e "The md5sum of the archive which contains all the above files is:" 
md5sum ${smrtcell_name}.${archive}.metadata.tar.gz
echo -e "##############################################################################\n"

echo -e "This is the links for the archive with the metadata and extra files\n"

s3cmd -c ~/.s3cfg_gentyane signurl s3://pacbio-data-provisoire/${mypath}/${smrtcell_name}.${archive}.metadata.tar.gz  $(echo "`date +%s` + 3600 * 24 * 14" | bc) |  sed "s/http/https/g"

rm ${smrtcell_name}.${archive}.metadata.tar.gz 


#Pour créer le md5 des bam :

for i in pacbio-data-provisoire/$mypath/hifi_reads/*.hifi_reads.[a-d]*.bam; do ln -sr $i .; done

echo -e "\n"

echo -e "##############################################################################"
echo -e "The md5sum of the BAM file(s) is:" 
#for i in *.hifi_reads.bc*.bam; do bash $HOME/compute_etag.sh $i 15; done
for i in *.hifi_reads.[a-d]*.bam; do md5sum $i; done
echo -e "##############################################################################\n"


#Créer liens pour les bam :
echo -e "This is the links for the BAM files\n"

for i in pacbio-data-provisoire/$mypath/hifi_reads/*.hifi_reads.[a-d]*.bam; do s3cmd -c ~/.s3cfg_gentyane signurl s3://$i  $(echo "`date +%s` + 3600 * 24 * 14" | bc); echo " "; done |  sed "s/http/https/g" 	


fusermount -u pacbio-data-provisoire

if [[ "$?" -eq "0" ]];then
	rm -rf $workingdir
fi
