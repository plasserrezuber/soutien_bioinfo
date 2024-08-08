#!/bin/bash

#SBATCH --job-name=hifiasm 
#SBATCH -o hifiasm-%A_%a.out
#SBATCH -n 192
#SBATCH --mem=512G
#SBATCH --partition=gdec
#SBATCH --time=3-00:00:00
#SBATCH --export=all

## SetUpFolders

workingdir=/storage/scratch/${USER}/${SLURM_JOBID}
TMPDIR=$workingdir/tmp
export TMPDIR 

mkdir -p $workingdir && cd $workingdir
pwd


echo -e "running on `hostname`\n"

## Set up the mount

module load python s3cmd goofys

export AWS_ACCESS_KEY_ID='b27b3d0a320941e3b4afab555faf0f2f'
export AWS_SECRET_ACCESS_KEY='f323f20b0dfd436295c44c0ded78dd99'

mkdir pacbio-data

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire ./pacbio-data


## symlinks for 3 ccs.bam

path1=r84087_20231025_074940/1_B01
path2=r84087_20231120_092219/1_B01
path3=r84087_20231123_081245/1_A01

ln -sr pacbio-data/$path1/hifi_reads/*hifi_reads.default.bam* ./
ln -sr pacbio-data/$path2/hifi_reads/*hifi_reads.default.bam* ./
ln -sr pacbio-data/$path3/hifi_reads/*hifi_reads.default.bam* ./

## convertir les bam en fasta.gz
ml smrttools/12.0.0.177059

bam2fasta -o ventricosa1 m84087_231025_082854_s2.hifi_reads.default.bam &
bam2fasta -o ventricosa2 m84087_231120_100037_s2.hifi_reads.default.bam &
bam2fasta -o ventricosa3 m84087_231123_082024_s1.hifi_reads.default.bam

## conda environment hifiasm
source /home/palasser/.bashrc
conda activate /home/palasser/apps/conda/envs/hifiasm_0.19

## assembly
## --write-ec: dump error corrected reads in FASTA
## --write-paf: read overlaps in PAF
## -l0 disable haplotig duplication purging
## --dual-scaf to output scaffoldings
## --n-hap: number of haplotype

hifiasm -o ventricosa_assembly -t192 -f39 --write-paf --write-ec --dual-scaf --n-hap 1 ventricosa*fasta.gz 


fusermount -u pacbio-data


# ### conversion gfa en fasta

#awk '/^S/{print ">"$2"\n"$3}' ventricosa_assembly.bp.p_ctg.gfa | fold > ventricosa_assembly.fa
grep 'S.*ptg' ventricosa_assembly.bp.p_ctg.gfa |cut -f2,3 |gawk '{print ">"$1"\n"$2}' > ae_ventricosa_assembly.fa

## 
ml samtools
samtools faidx ae_ventricosa_assembly.fa

# ## stat sur assemblage
ml genometools

gt seqstat ae_ventricosa_assembly.fa > ae_ventricosa_seqstat.txt

## copy results
cd $HOME
mv $workingdir/ventricosa_assembly*  /home/palasser/projects/soutien_bioinfo/SOURDILLE/results/   #&& rm -rf $workingdir

