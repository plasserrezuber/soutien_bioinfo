#!/bin/bash

#SBATCH --job-name=archivage
#SBATCH -o archivage-%A_%a_task.out
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --partition=normal
#SBATCH --export=all


if [ "$#" -ne 1 ]; then
        echo -e "Need 2 arguments :\n"
        echo -e "-> The path of the files to archive without the   /   at the end : eg. /r64071_20200805_142732/2_B01 \n"
	exit
fi

echo -e "running on `hostname`\n"

#SetUpFolders


workingdir='/storage/scratch/'$SLURM_JOBID

mkdir $workingdir && cd $workingdir && mkdir subreads

echo -e "The workingdir's path is :" $(pwd) "\n"


#LoadModules

module load python s3cmd goofys

mypath=${1?Usage: $0 <mypath> <smrtcell_name>}

#Set up the mount

#Set up the name of the archive

#echo -e "Ok, the archive ${archive}.tar is gonna be created with the files in $mypath ; and then that archive will be put in the project $project in OSCAR\n"


### To check if the scraps files are already deleted or not ###

echo "########################################################"
echo "Checking if the scraps files are deleted or not in OSCAR"
echo -e "########################################################\n"

if [[ $(s3cmd -c ~/.s3cfg_gentyane del s3://pacbio-data-provisoire/$mypath/*scraps*) == *delete* ]]; then
        echo -e "The scraps files are now deleted before making the archive\n";
else
        echo -e "The scraps files were already deleted\n";
fi

### In order to create an archive, be sure to make hidden files visible ###

echo -e "#############################################"
echo -e "Just moving the hidden files in visible files"
echo -e "#############################################\n"

if [[ $(s3cmd -c ~/.s3cfg_gentyane get s3://pacbio-data-provisoire/$mypath/.*) == download* ]]; then
        for file in .[^.]*;
        do
                if [[ -f "$file" ]]; then
                mv "${file}" "${file#.}";
                fi;
        done
        echo -e "Removing hidden files in s3://pacbio-data/$mypath/\n"
        s3cmd -c ~/.s3cfg_gentyane del s3://pacbio-data-provisoire/$mypath/.*
        echo -e "\n";
        echo -e "Adding visible files into s3://pacbio-data/$mypath/\n"
        s3cmd -c ~/.s3cfg_gentyane put *metadata.xml s3://pacbio-data-provisoire/$mypath/
        echo -e "\n";
        rm *metadata.xml;
else
        echo -e "There are no hidden files s3://pacbio-data/$mypath/\n"
fi


#Set up the mount

mkdir pacbio-data

export AWS_ACCESS_KEY_ID='b27b3d0a320941e3b4afab555faf0f2f'
export AWS_SECRET_ACCESS_KEY='f323f20b0dfd436295c44c0ded78dd99'

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire ./pacbio-data


if compgen -G "pacbio-data/$mypath/*subreadset.xml" > /dev/null; then
        archive=$(fgrep 'PacBio.DataSet.SubreadSet' pacbio-data/$mypath/*subreadset.xml |cut -d '=' -f 4| cut -d '"' -f 2)
        archive="${archive// /_}"
	archive="${archive//-Cell[1-9]/}" 
	smrtcell_name=$(fgrep 'ResourceId' pacbio-data/$mypath/*subreadset.xml  | cut -d '=' -f 4 |grep '.subreads.bam.pbi' |cut -d '"' -f 2 |cut -d '.' -f1) 
        #project=$(fgrep 'PacBio.DataSet.SubreadSet' pacbio-data/$mypath/*subreadset.xml |cut -d '=' -f 4| cut -d '"' -f 2| cut -d '-' -f 1)
fi

echo -e "\n"

### Change the subreads.bam name with its real name ###
#s3cmd -c ~/.s3cfg_gentyane mv s3://pacbio-data/$mypath/${smrtcell_name}.subreads.bam s3://pacbio-data/$mypath/$archive.subreads.bam

echo -e "Creating and adding the md5sum file to the bucket"

### compute the md5 du fichier $archive.subreads.bam###
s3cmd -c ~/.s3cfg_gentyane ls --list-md5 s3://pacbio-data-provisoire/$mypath/${smrtcell_name}.subreads.bam > ${smrtcell_name}.subreads.md5 

s3cmd -c ~/.s3cfg_gentyane put *.subreads.md5 s3://pacbio-data-provisoire/$mypath/${smrtcell_name}.subreads.md5

echo -e "Done\n"

### Create the symlinks and the archive ###

fusermount -u pacbio-data

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire ./pacbio-data

echo -e "##########################################################################################"
echo -e "Creation of the symlinks for the files present in pacbio-data/$mypath"
echo -e "##########################################################################################\n"

for i in pacbio-data/$mypath/*;
do if [[ -f $i ]]; then
   ln -sr ${i} ./subreads/;
   fi;
done

if [[ "$?" -eq "0" ]];then
        echo -e "The symlinks are now created\n"
fi

cd subreads


echo -e "#########################################"
echo -e "Creation of the archive ${smrtcell_name}.metadata.tar.gz"
echo -e "#########################################\n"

echo -e "The archive ${smrtcell_name}.${archive}.metadata.tar.gz will contain theses files :\n"

tar -czvhf ${smrtcell_name}.${archive}.metadata.tar.gz * --exclude='*.subreads.bam'

echo -e "\n"

cd ..


echo -e "##############################################################################"
echo -e "Transfer of ${smrtcell_name}.${archive}.metadata.tar.gz on OSCAR in the bucket pacbio_data-provisoire/$mypath"
echo -e "##############################################################################\n"

fusermount -u pacbio-data

if [[ "$?" -eq "0" ]];then
        s3cmd -c ~/.s3cfg_gentyane put subreads/${smrtcell_name}.${archive}.metadata.tar.gz s3://pacbio-data-provisoire/${mypath}/${smrtcell_name}.${archive}.metadata.tar.gz
        echo -e "\n"
        cd $HOME
        rm -rf $workingdir
fi

echo -e "The archive ${smrtcell_name}.${archive}.metadata.tar.gz is now on OSCAR\n"

echo -e "And this is the links for the subreads\n" 

s3cmd -c ~/.s3cfg_gentyane signurl s3://pacbio-data-provisoire/${mypath}/${smrtcell_name}.subreads.bam $(echo "`date +%s` + 3600 * 24 * 14" | bc) | sed "s/http/https/g"

echo -e "\n" 

echo -e "And this is the links for the metadata\n"

s3cmd -c ~/.s3cfg_gentyane signurl s3://pacbio-data-provisoire/${mypath}/${smrtcell_name}.${archive}.metadata.tar.gz $(echo "`date +%s` + 3600 * 24 * 14" | bc) | sed "s/http/https/g" 

echo -e "\n" 

sacct -o reqmem,maxrss,maxvmsize,averss,elapsed,alloccpus -j $SLURM_JOB_ID 
