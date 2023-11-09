#!/bin/bash

#SBATCH --job-name=ccs 
#SBATCH -o ccs-%A_%a_task.out
#SBATCH -n 128
#SBATCH --mem=64G
#SBATCH --partition=gdec
#SBATCH --time=3-00:00:00
#SBATCH --export=all

#SetUpFolders

workingdir='/storage/scratch/'$SLURM_JOBID
TMPDIR=$workingdir/tmp
export TMPDIR 

mkdir $workingdir && cd $workingdir && mkdir subreads
pwd

mypath=${1?Usage: $0 <mypath> <mydir>}
mydir=${2?Usage: $0 <mypath> <mydir>}

mkdir -p $mydir

echo -e "Ok, the CCS are gonna be created with the subreads file in pacbio-data/$mypath\n" 

#LoadModules

module load python s3cmd goofys

echo -e "running on `hostname`\n"

#Set up the mount

export AWS_ACCESS_KEY_ID='b27b3d0a320941e3b4afab555faf0f2f'
export AWS_SECRET_ACCESS_KEY='f323f20b0dfd436295c44c0ded78dd99'

mkdir pacbio-data

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire ./pacbio-data

### To check if the scraps files are already deleted or not ###

echo "###########################################################"
echo "Checking if the scraps files are deleted or not in OSCAR"
echo -e "########################################################\n"

if [[ $(s3cmd -c ~/.s3cfg_gentyane del s3://pacbio-data-provisoire/$mypath/*scraps*) == *delete* ]]; then
        echo -e "The scraps files are now deleted before making the archive\n";
else
        echo -e "The scraps files were already deleted\n";
fi

### Create the symlinks and the archive ###

echo -e "##########################################################################################"
echo -e "Creation of the symlinks for the subreads files present in pacbio-data/$mypath"
echo -e "##########################################################################################\n"

for i in pacbio-data/$mypath/*subreads.bam*;
do if [[ -f $i ]]; then
   ln -sr ${i} ./subreads/;
   fi;
done

if [[ "$?" -eq "0" ]];then
        echo -e "The symlinks are now created\n"
fi

echo -e "##########################################################################################"
echo -e "Generation of the CCS from the subreads file in pacbio-data/$mypath"
echo -e "##########################################################################################\n"

cd subreads

module load smrttools/12.0.0.177059

echo "start ccs at `date`"

echo "ccs -j 128 --min-passes 3 --min-rq 0.99 *.subreads.bam ccs.bam"

time ccs -j 128  --min-rq 0.99 --min-passes 3  *.subreads.bam ccs.bam

echo -e "\n"

echo "stop ccs at `date`"

echo -e "\n"

echo -e "generate runqc plots of the ccs file\n"

dataset create ccs.consensusreadset.xml ccs.bam
runqc-reports ccs.consensusreadset.xml 

rm *thumb*

echo -e "\n"

module load samtools genometools

samtools fasta ccs.bam > ccs.fa
gt seqstat ccs.fa > ccs_stats.txt
sed -i 's/contigs/CCS/g; s/contig/CCS/g' ccs_stats.txt
rm ccs.fa

cd ..

fusermount -u pacbio-data


if [[ "$?" -eq "0" ]];then
        mv $workingdir/subreads/ccs* $mydir/ 
        cd $HOME
        rm -rf $workingdir
fi

echo -e "\n"

echo -e "The CCS are now in $mydir" 

echo -e "\n" 

sacct -o reqmem,maxrss,maxvmsize,averss,elapsed,alloccpus -j $SLURM_JOB_ID

