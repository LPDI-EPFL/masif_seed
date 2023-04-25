#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 4096
#SBATCH --time 01:30:00
#SBATCH --array=1-363
#SBATCH --output=exelogs/out/splitfile.%A_%a.out
#SBATCH --error=exelogs/err/splitfile.%A_%a.err

# To be modified in this file
# 1) The path ('in_dir') to the database of scaffolds for grafting
# 2) The number of arrays for paralellization in regards to the number of splitfiles that composed the scaffold database
# 3) The location of the Rosetta binary files

in_dir=/path/to/scaffold/database
ROSETTA_BIN=./rosetta_07_2020/main/source/bin/rosetta_scripts.linuxiccrelease

MYID=${SLURM_ARRAY_TASK_ID}
UPPERID=$1
HOTSPOTS=$2

mkdir -p results/seed_$UPPERID/$MYID
$ROSETTA_BIN -overwrite \
         -ex1 -ex2 -ignore_zero_occupancy false \
         -parser:protocol grafting_varhotspots.xml\
         -l $in_dir/splitfile_$MYID.split \
         -ignore_unrecognized_res -nstruct 1 \
         -out:path:all results/seed_$UPPERID/$MYID \
         -parser:script_vars r=$HOTSPOTS s=$UPPERID c=$UPPERID \
         -no_nstruct_label \
         -out:suffix _S$UPPERID

echo "CASTOR: RUN FINISHED"
