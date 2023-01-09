#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32768
#SBATCH --time 01:30:00
#SBATCH --array=1-363
#SBATCH --output=exelogs/out/splitfile.%A_%a.out
#SBATCH --error=exelogs/err/splitfile.%A_%a.err

in_dir=/work/lpdi/users/anthony/database/small_prots
ROSETTA_BIN=/work/lpdi/bin/rosetta_07_2020/main/source/bin/rosetta_scripts.linuxiccrelease

MYID=${SLURM_ARRAY_TASK_ID}
UPPERID=$1
HOTSPOTS=$2
NFRAGS=$3

if [ $NFRAGS = 1 ]; then
    DELTA1=0:0
    DELTA2=0:0
elif [ $NFRAGS = 2 ]; then
    DELTA1=0:0,0:0
    DELTA2=0:0,0:0
elif [ $NFRAGS = 3 ]; then
    DELTA1=0:0,0:0,0:0
    DELTA2=0:0,0:0,0:0
else
    echo "Error with the number of fragments"
fi

mkdir -p results/seed_$UPPERID/$MYID
$ROSETTA_BIN -overwrite \
         -ex1 -ex2 -ignore_zero_occupancy false \
         -parser:protocol grafting_varhotspots2.xml\
         -l $in_dir/splitfile_$MYID.split \
         -ignore_unrecognized_res -nstruct 1 \
         -out:path:all results/seed_$UPPERID/$MYID \
         -parser:script_vars h=$HOTSPOTS s=$UPPERID c=$UPPERID f=$DELTA1 r=$DELTA2\
         -no_nstruct_label \
         -out:suffix _S$UPPERID
#         -parser:protocol minimal_motif_graft_and_design.xml\
echo "CASTOR: RUN FINISHED"
