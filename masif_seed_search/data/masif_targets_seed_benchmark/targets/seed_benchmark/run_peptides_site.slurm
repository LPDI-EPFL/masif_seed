#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32000
#SBATCH --time 1:30:00
#SBATCH --array=100-130

ppi_pair_id=$(cat lists/split/x${SLURM_ARRAY_TASK_ID})
time ./run.sh $ppi_pair_id
