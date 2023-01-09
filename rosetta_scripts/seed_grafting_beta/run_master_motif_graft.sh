#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4096
#SBATCH --time 00:10:00
#SBATCH --array=108,109,112,114,115,116,118,119,120,122,123,124,127,129,132,135,136,138,139,140,143,147,149,150
#SBATCH --output=exelogs/out/splitfile.%A_%a.out
#SBATCH --error=exelogs/err/splitfile.%A_%a.err

MYID=${SLURM_ARRAY_TASK_ID}
COMPLEX=$(sed -n "$SLURM_ARRAY_TASK_ID"p list.txt)
python3 /work/lpdi/users/anthony/seed_grafting/tools/extract_complex_multifrag.py $COMPLEX $MYID
python3 /work/lpdi/users/anthony/seed_grafting/tools/crop_seed.py $MYID
HOTSPOTS=$(python3 /work/lpdi/users/anthony/seed_grafting/tools/find_hotspots_multifrag.py $MYID)
NFRAGS=$(python3 /work/lpdi/users/anthony/seed_grafting/tools/find_nfrag.py $MYID)
echo "Hotspots: "$HOTSPOTS
echo "Number of fragments: "$NFRAGS
mkdir -p results/seed_$MYID

sbatch run_motif_graft_varhotspots.slurm $MYID $HOTSPOTS $NFRAGS

echo "CASTOR: RUN FINISHED"
