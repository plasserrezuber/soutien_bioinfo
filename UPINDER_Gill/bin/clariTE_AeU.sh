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
OUTPUT='/home/palasser/projects/soutien_bioinfo/UPINDER_Gill/results/annotTE_AeU/'$1

mkdir -p -m 750 $OUTPUT
cd $OUTPUT

files=($(find $OUTPUT/$1*.fa.out.xm))

i=${files[$SLURM_ARRAY_TASK_ID]}

i_fasta=$(cut -d'.' -f1,2 <(echo $i))

# mise en forme
sed -i "s/#Unspecified/#Unknown/" $i

# ##### clariTE
$DATABANK/bin/clariTE.pl -dir $OUTPUT/ \
-LTR $DATABANK/databank/CLARIwheat.LTR_position \
-classi $DATABANK/databank/CLARIwheat.classification \
-fasta $i_fasta \
$i
