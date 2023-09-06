#!/bin/bash

###### embl2gff: takes all chunks embl and produces one merged .gff per chrom
## NB: script embl2gff.generic.pl gives an ID per TE per chrom (not per chunk) from zero if a list of ordered files for all chunks is given for each chrom
## That is why the $embl_files variable is created before

embl_files=$(\ls -1 ../results/${1}/${1}*.embl |sort -t ':' -k2,2n |tr -s '\n' ' ')
## !!!!! revoir option "renan"
bin/embl2gff.generic.pl -renan -note -l 10 $embl_files > ../results/${1}_clariTE.gff

###### gff to gff3: new coordinates calculation with gawk command (from chunks relative to chrom relative)
grep -v $'\t''region' ../results/${1}_clariTE.gff |sed 's/embl2gff.generic.pl-2.1/clariTE/' \
|gawk -v FS='\t' -v OFS='\t' '/^chr/ {match($1, /:[0-9]+-/); $4=$4+substr($0, RSTART+1, RLENGTH-2); $5=$5+substr($0, RSTART+1, RLENGTH-2); print}' \
|gt gff3 -sort -tidy -retainids 1> ../results/${1}_clariTE_tmp.gff3  2> ../results/${1}_clariTE_tmp_gt.log

### variable $endchrom to give the right sequence-region coordinate in gff3
endchrom=$(grep ${1} ${2}.fai |cut -f2)

# grep -v command keeps only one "sequence-region" line starting with zero coordinate
# gawk and sed commands to format in a more friendly way
grep -v -P '##sequence-region *'$1':[1-9]' ../results/${1}_clariTE_tmp.gff3 \
|gawk -v LG=$1 -v end=$endchrom 'BEGIN{FS='\t';OFS='\t'} { if ($0~'##sequence-region') $0='##sequence-region\t'LG'\t1\t'end; print }' \
|sed -E 's/'$1':[0-9]*-[0-9]*/'$1'/' |sed -E 's/Compo:.* (Family)/\1/' | sed -E 's/Post:.* (Status)/\1/' > ../results/${1}_clariTE.gff3
