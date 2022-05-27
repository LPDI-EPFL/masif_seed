#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --partition serial
#SBATCH --cpus-per-task 1
#SBATCH --mem 4096
#SBATCH --time 0:10:00
#SBATCH --array=1-2545
echo STARTING AT `date`
echo ${SLURM_ARRAY_TASK_ID}
echo FINISHED at `date`

PDBID_CHAIN=$1
masif_seed_search_root=$(git rev-parse --show-toplevel)
masif_seed_search_scripts=$masif_seed_search_root/scripts
export PYTHONPATH=$PYTHONPATH:$masif_seed_search_scripts:`pwd`
echo $masif_seed_search_scripts
python3 -W ignore $masif_seed_search_scripts/align_all_to_all_fixed_length.py $PDBID_CHAIN ${SLURM_ARRAY_TASK_ID}


