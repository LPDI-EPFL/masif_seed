#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16384
#SBATCH --time 23:30:00
#SBATCH --array=1-31
echo STARTING AT `date`
echo ${SLURM_ARRAY_TASK_ID}

i=1
while read p; do
    if [ $(( i % 32 )) == ${SLURM_ARRAY_TASK_ID} ]; then
        ./run_one.sh $p
    fi
    i=$((i+1))
done < ../benchmark_lists/merged_receptor_helix_list.txt


echo FINISHED at `date`
