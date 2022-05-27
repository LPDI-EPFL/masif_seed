#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16384
#SBATCH --time 23:0:00
#SBATCH --array=1-100

# Run ZDock without blocked residues, only count the time of this step.
# Now run the target against each decoy (which will include the real binder).
i=1
PDBID=$1
CHAIN1=$2
while read decoy_ppi_pair
do
    if [ $(( i % 101 )) == ${SLURM_ARRAY_TASK_ID} ]; then
	DECOY_PDBID=$(echo $decoy_ppi_pair | cut -d"_" -f1) 
	DECOY_CHAIN=$(echo $decoy_ppi_pair | cut -d"_" -f2) 
	decoy_filename=../../01-zdock_marked/$DECOY_PDBID\_$DECOY_CHAIN\_m.pdb
		
	runname=$PDBID\_$CHAIN1\_$DECOY_PDBID\_$DECOY_CHAIN
	# Run ZDock with blocked residues, only count the time of this step.
	/usr/bin/time -f 'user %U' -o $runname\_cpu_seconds.txt ./zdock -S 100 -o zdock_$runname\.out -R $PDBID\_$CHAIN1\_m\_bl.pdb -L $decoy_filename
	# Run ZDock without blocked residues, only count the time of this step.
#        /usr/bin/time -f 'user %U' -o $runname\_cpu_seconds.txt ./zdock -S 100 -o zdock_$runname\.out -R $PDBID\_$CHAIN1\_m.pdb -L $decoy_filename
    fi
    i=$((i+1))
done < input_list.txt


