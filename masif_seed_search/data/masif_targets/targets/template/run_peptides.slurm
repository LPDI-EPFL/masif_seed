#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16384
#SBATCH --time 1:30:00
#SBATCH --job-name=""
#SBATCH --array=1-1

if [[ $# -eq 0 ]] ; then
    echo 'You must pass as an argument the name of the pdb and the chain (e.g. 3R2X_AB)'
    exit 0
fi

source /work/upcorreia/bin/load_masif_environment.sh

masif_seed_search_root=$(git rev-parse --show-toplevel)
masif_root=$masif_seed_search_root/../
masif_target_root=$masif_seed_search_root/data/masif_targets/
export masif_db_root=../../../../../masif/
masif_source=$masif_root/source/
masif_data=$masif_root/data/
export masif_root
export masif_target_root
export PYTHONPATH=$PYTHONPATH:$masif_source:`pwd`
python -W ignore $masif_seed_search_root/source/masif_seed_search_nn.py params_peptides $1 lists/split_peptides/${SLURM_ARRAY_TASK_ID}
