#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4096
#SBATCH --time 00:10:00
#SBATCH --array=1-25
#SBATCH --output=exelogs/out/splitfile.%A_%a.out
#SBATCH --error=exelogs/err/splitfile.%A_%a.err

# To be modified in this file : The number of arrays ('#SBATCH --array=...') corresponding to the number of seeds to be grafted

COMPLEX=$(sed -n "$SLURM_ARRAY_TASK_ID"p list.txt)
CHAINS=($(python3 ../tools/find_chains.py $COMPLEX))
HOTSPOTS=$(python3 ../tools/find_hotspots.py $COMPLEX ${CHAINS[0]} ${CHAINS[1]})
MYID=${SLURM_ARRAY_TASK_ID}

python3 ../tools/extract_complex.py $COMPLEX $MYID
echo $HOTSPOTS
mkdir -p results/seed_$MYID

sbatch run_motif_graft_varhotspots.slurm $MYID $HOTSPOTS

echo "CASTOR: RUN FINISHED"
