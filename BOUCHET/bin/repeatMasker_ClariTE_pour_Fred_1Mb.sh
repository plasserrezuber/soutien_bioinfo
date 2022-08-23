#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --partition=fast
#SBATCH --qos=fast
#SBATCH --export=ALL
#SBATCH -J arLrTE

module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 RepeatMasker bioperl/1.7.0_rc5

INPUT='/home/frchoule/data/projects/rustwatch/markers/contextSeq_round2/Lr14a_MT123593/alignment_1Mb_region'
OUTPUT='/home/palasser/soutien_bioinfo/BOUCHET/annot_TEs_arinaLrFor'
DATABANK='/storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0'

mkdir -p -m 750 $OUTPUT
cd $OUTPUT

cp $INPUT/arinaLrFor_Lr14a_1MbRegion.fa $OUTPUT/
###################################################################################
##### RepeatMasker
# -xsmall returns repetitive regions in lowercase (rest capitals) rather than masked
# -norna interested in small RNA genes (mostly tRNAs and snRNAs), you should use the -norna option that leaves these sequences unmasked, while still masking SINEs.
# -nolow         does not mask low complexity DNA or simple repeats
# -cutoff [number] sets cutoff score for masking repeats when using -lib (default cutoff 225)

RepeatMasker -e crossmatch -lib $DATABANK/databank/CLARIwheat.fasta -xsmall -nolow -xm -pa $SLURM_CPUS_PER_TASK -q $OUTPUT/arinaLrFor_Lr14a_1MbRegion.fa


sed -i "s/#Unspecified/#Unknown/" $OUTPUT/arinaLrFor_Lr14a_1MbRegion*.out.xm

# ##### clariTE
$DATABANK/bin/clariTE.pl -dir $OUTPUT/ \
-LTR $DATABANK/databank/CLARIwheat.LTR_position \
-classi $DATABANK/databank/CLARIwheat.classification \
-fasta $OUTPUT/arinaLrFor_Lr14a_1MbRegion.fa \
$OUTPUT/arinaLrFor_Lr14a_1MbRegion.fa.out.xm