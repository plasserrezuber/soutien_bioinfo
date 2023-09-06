#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH --partition=gdec
#SBATCH --export=ALL

module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker/4.0.5

OUTPUT='/home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7/'$1
SCRATCHDIR='/storage/scratch/'$USER'/Lo7/'$1

mkdir -p -m 750 $OUTPUT
mkdir -p -m 750 $SCRATCHDIR
cd $SCRATCHDIR

files=($(find $OUTPUT/$1*.fa))
i=${files[$SLURM_ARRAY_TASK_ID]}

prefix=$(basename $i)

cp $i $SCRATCHDIR/$prefix

###################################################################################
##### RepeatMasker
# -pa parallel
# -xm     creates an additional output file in cross_match format (for parsing)
# -q             Quick search; 5-10% less sensitive, 3-4 times faster than default
# -xsmall returns repetitive regions in lowercase (rest capitals) rather than masked
# -norna interested in small RNA genes (mostly tRNAs and snRNAs), you should use the -norna option that leaves these sequences unmasked, while still masking SINEs.
# -nolow         does not mask low complexity DNA or simple repeats
# -cutoff [number] sets cutoff score for masking repeats when using -lib (default cutoff 225)

RepeatMasker -e crossmatch -lib /storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0/databank/CLARIwheat.fasta \
-xsmall -nolow -xm -pa $SLURM_CPUS_PER_TASK -q $SCRATCHDIR/$prefix

mv $SCRATCHDIR/${prefix}.out.xm $OUTPUT/ && rm $SCRATCHDIR/${prefix}*