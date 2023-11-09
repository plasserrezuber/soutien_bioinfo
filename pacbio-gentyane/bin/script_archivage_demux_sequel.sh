#!/bin/bash

#SBATCH --job-name=demux
#SBATCH -o archivage-%A_task.out
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --partition=gdec
#SBATCH --time=03:00:00
#SBATCH --export=all


### ce script prend 1 argument

### variables
myproject=${1?Usage: $0 <myproject>}


mypath='/home/palasser/projects/soutien_bioinfo/pacbio-gentyane/results'


### dossier de travail
mkdir -p ${mypath}/${myproject}/demultiplexed_bam

cd ${mypath}/${myproject}

echo -e "The workingdir's path is :" $(pwd) "\n"

##### demultiplexage   =============> CHANGER NOMBRE CPU SI LIMA 

### le script prendrait alors 2 arguments
#barcode_fasta=${2?Usage: $0 <mypath> <barcode_fasta>}

#module load smrttools/12.0.0.177059
#lima -j 4 --ccs --peek-guess --split-bam-named ccs.bam ${barcode_fasta} ./demultiplexed_bam/demux.bam

##### archivage

tar --exclude='ccs.bam' --exclude='*demu*' --exclude='slurm*' --exclude='*.metadata.tar.gz' -czvf ${myproject}.metadata.tar.gz ./

echo -e "\n"

tar -cvf demultiplexed_bam.tar ./demultiplexed_bam

echo -e "\n"

#LoadModules
module load python s3cmd 

#### transfert sur Oscar

s3cmd -c ~/.s3cfg_gentyane put ./ccs.bam s3://gentyane-share/${myproject}/ccs.bam
s3cmd -c ~/.s3cfg_gentyane put ./demultiplexed_bam.tar s3://gentyane-share/${myproject}/demultiplexed_bam.tar
s3cmd -c ~/.s3cfg_gentyane put ./${myproject}.metadata.tar.gz s3://gentyane-share/${myproject}/${myproject}.metadata.tar.gz


##### liens telechargement 

echo -e "##############################################################################"
echo -e "The md5sum of the metadata archive is:" 
md5sum ${myproject}.metadata.tar.gz
echo -e "##############################################################################\n"

echo -e "This is the link for the archive with the metadata and extra files\n"

s3cmd -c ~/.s3cfg_gentyane signurl s3://gentyane-share/${myproject}/${myproject}.metadata.tar.gz $(echo "`date +%s` + 3600 * 24 * 14" | bc) |  sed "s/http/https/g"

echo -e "\n"

echo -e "##############################################################################"
echo -e "The md5sum of the ccs.bam file is:" 
md5sum ccs.bam
echo -e "##############################################################################\n"

echo -e "This is the link for the ccs.bam\n"

s3cmd -c ~/.s3cfg_gentyane signurl s3://gentyane-share/${myproject}/ccs.bam $(echo "`date +%s` + 3600 * 24 * 14" | bc) |  sed "s/http/https/g"

echo -e "\n"

echo -e "##############################################################################"
echo -e "The md5sum of the demultiplexed_bam archive is:" 
md5sum demultiplexed_bam.tar
echo -e "##############################################################################\n"

echo -e "This is the link for the demultiplexed bam\n"

s3cmd -c ~/.s3cfg_gentyane signurl s3://gentyane-share/${myproject}/demultiplexed_bam.tar $(echo "`date +%s` + 3600 * 24 * 14" | bc) |  sed "s/http/https/g"

echo -e "\n"