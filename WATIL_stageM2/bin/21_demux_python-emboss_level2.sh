#!/bin/bash
#SBATCH --job-name=demux2
#SBATCH -p fast
#SBATCH --qos=fast
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --array=0-93

n=1
for files in demux.bc*fasta; do array[n++]+=$files; done

module load python/3.4.5

echo "start demultiplexing at `date`"

for i in $(echo "${array[$SLURM_ARRAY_TASK_ID]}");
do
python /home/newatil/bingit/Demultiplexing_v2.py ${i} /home/newatil/data/E779-Primers.fasta
done

echo "end demultiplexing at `date`"


############################################################################
##EXECUTE DANS TERMINAL
#GENES=$(fgrep '>' /home/newatil/data/E779-Primers.fasta |sed 's/>Glu_//' |cut -d'-' -f1 |sort |uniq)

#for g in $GENES; do len=$(cat /home/newatil/python-emboss/demux*${g}*fasta.fai |cut -f2 |sort -n |uniq -c |sed 's/ */\t/' |tr ' ' '\t' |cut -f2,3 |sort -k1,1n |tail -n1 |cut -f2); for f in $(find /home/newatil/python-emboss/demux.*_Glu_${g}_test.fasta); do f2=$(basename $f |cut -d'.' -f1,2); awk -v l=$len -v OFS='\t' '{ if ($2>l-1000 && $2<l+1000) print $1,"0",$2,$1}' ${f2}.fasta.fai > ${f2}_selection.bed; bedtools getfasta -fi $f -bed ${f2}_selection.bed -fo ${f2}_selection.fasta; done; done
############################################################################