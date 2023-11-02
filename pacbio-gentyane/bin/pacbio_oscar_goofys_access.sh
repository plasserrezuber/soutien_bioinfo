#!/bin/bash

module load python goofys s3cmd


#connexion s3 avec goofys= montage du bucket contenant data a utiliser
#s3cmd ls s3://wheatomics-raw/assemblies/Triticum_aestivum_RENAN_v1.fasta

export AWS_ACCESS_KEY_ID=b27b3d0a320941e3b4afab555faf0f2f
export AWS_SECRET_ACCESS_KEY=f323f20b0dfd436295c44c0ded78dd99

mkdir $HOME/s3-pacbio-gentyane

goofys --stat-cache-ttl 3600s -uid $UID -gid $(id -g) --dir-mode=0500 --file-mode=0400 --cheap --endpoint https://s3.mesocentre.uca.fr pacbio-data-provisoire $HOME/s3-pacbio-gentyane


####################################################################################
#couper la connection= defaire le montage du bucket du CEPH
#fusermount -u $HOME/s3-pacbio-gentyane
