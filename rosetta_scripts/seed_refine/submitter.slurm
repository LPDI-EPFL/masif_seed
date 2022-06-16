#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8192
#SBATCH --time 1:00:00
#SBATCH --array=1-2524
#SBATCH --output=./exelogs/out/splitfile.%A_%a.out
#SBATCH --error=./exelogs/err/splitfile.%A_%a.err

# To be modified in this file: 
# 1) Path to the Rosetta binary files
# 2) Number of arrays (= number of seeds) in "#SBATCH --array=..."
# To be provided in the working repository :
# 1) A list of the seed and target filename
# 2) The input PDB of the seeds in the "in" directory

ROSETTA_BIN=/rosetta_07_2020/main/source/bin/rosetta_scripts.linuxiccrelease

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p list.txt)
MYID=$SLURM_ARRAY_TASK_ID
SEED_PATH=$(echo $LINE | awk '{print $1;}')
SEED_NAME=$(basename $SEED_PATH | cut -d "." -f1)
IN_FILE=./relax/$SEED_NAME*.pdb
echo $LINE > ./tmp/tmp_list_$MYID.txt | $LINE

$ROSETTA_BIN -parser:protocol relax.xml \
    -list ./tmp/tmp_list_$MYID.txt \
    -no_nstruct_label \
    -out:path:all ./relax \
    -overwrite \
    -renumber_pdb \
    -ignore_unrecognized_res

$ROSETTA_BIN -parser:protocol interface_design.xml \
    -s $IN_FILE \
    -no_nstruct_label \
    -out:path:all ./out \
    -overwrite \
    -renumber_pdb \
    -ignore_unrecognized_res \

rm ./tmp/tmp_list_$MYID.txt
