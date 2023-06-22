#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --partition=gdec
#SBATCH --export=ALL

module load gcc/4.8.4 triannotTools/1.2 gdecTools/1.1 bioperl/1.7.0_rc5

##### DEFINITION VARIABLES
DATABANK='/storage/groups/gdec/triannot/tools/ClariTE/clariTE_1.0'
OUTPUT='/home/palasser/projects/soutien_bioinfo/BOUGUENNEC/results/annotTE_Lo7/'$1

mkdir -p -m 750 $OUTPUT
cd $OUTPUT

files=($(find $OUTPUT/$1*.fa))
i=${files[$SLURM_ARRAY_TASK_ID]}

# mise en forme
sed -i "s/#Unspecified/#Unknown/" $i.out.xm

# ##### clariTE
$DATABANK/bin/clariTE.pl -dir $OUTPUT/ \
-LTR $DATABANK/databank/CLARIwheat.LTR_position \
-classi $DATABANK/databank/CLARIwheat.classification \
-fasta $i \
$i.out.xm
