#!/bin/bash

#SBATCH --job-name=demuxrev
#SBATCH -o demuxrev-%A_%a_task.out
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --partition=gdec
#SBATCH --time=3-00:00:00
#SBATCH --export=all

#SetUpFolders

workingdir='/storage/scratch/'$SLURM_JOBID
TMPDIR=$workingdir/tmp
export TMPDIR 

mkdir $workingdir && cd $workingdir
pwd

#LoadModules

module load python s3cmd goofys

echo -e "running on `hostname`\n"

#Set up the mount

export AWS_ACCESS_KEY_ID='b27b3d0a320941e3b4afab555faf0f2f'
export AWS_SECRET_ACCESS_KEY='f323f20b0dfd436295c44c0ded78dd99'


mkdir pacbio-data

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire ./pacbio-data

mypath='r84087_20231120_092219/1_D01'
mybam='m84087_231120_110249_s4.hifi_reads.unassigned.bam'
cp /home/palasser/projects/soutien_bioinfo/pacbio-gentyane/data/E1212_barcodes.fasta ./ 

ln -sr pacbio-data/${mypath}/hifi_reads/${mybam} ./

module load smrttools/12.0.0.177059

lima -j 16 --ccs --peek-guess --split-bam-named ${mybam} E1212_barcodes.fasta m84087_231120_110249_s4.hifi_reads_demux.bam

s3cmd -c ~/.s3cfg_gentyane put *demux* s3://pacbio-data-provisoire/${mypath}/hifi_reads/

fusermount -u pacbio-data


rm -rf $workingdir
