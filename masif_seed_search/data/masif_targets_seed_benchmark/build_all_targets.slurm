#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16384
#SBATCH --time 2:0:0
#SBATCH --array=2-32
echo STARTING AT `date`
echo ${SLURM_ARRAY_TASK_ID}
echo FINISHED at `date`
source /work/upcorreia/bin/load_masif_environment.sh
i=1
while read p; do
    if [ $(( i % 1001 )) == ${SLURM_ARRAY_TASK_ID} ]; then
        ./run_target_protocol.sh $p
    fi
    i=$((i+1))
done < lists/nonredundant_helices_extended_receptoronly.txt

