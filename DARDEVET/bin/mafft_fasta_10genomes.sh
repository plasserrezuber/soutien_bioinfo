# !/bin/bash
# #SBATCH --time=1:00:00 #1h
# SBATCH --nodes=1 # Un noeud par tache
# SBATCH --ntasks=1
# SBATCH --mem=8G
# SBATCH --cpus-per-task=8  # Nb of threads we want to run on (il y a 32 CPU/noeud)
# SBATCH --job-name=mafft
# SBATCH --partition=smp

# ce script peut etre lance sur le frontal car demande peut de ressources

module restore
module swap gcc/4.8.4 gcc/8.1.0
module load openmpi/3.0.0 MAFFT/7.427

for file in $(find $HOME/results/bedtools/tenGenomes/bygene -not -empty -name '*.fna')
do
  mafft --quiet --globalpair --maxiterate 1000 $file > $HOME/results/mafft/$(basename $file .fna)_clustal.fna || echo "mafft failed for $file"
done
