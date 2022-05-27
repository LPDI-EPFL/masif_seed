#!/bin/sh

#for ppi_pair in 3F74_A_B 3HRD_E_H 1ERN_A_B 2B3Z_C_D 4KGG_C_A 2FE8_A_C 1SOT_A_C 2LBU_E_D 3PGA_1_4 3QWN_I_J 3TND_B_D 1XDT_T_R
NUM_DECOYS=2000
helicesdir=/home/gainza/lpdi_fs/masif_design_paper/repository/masif_seed/masif/data/masif_helices_benchmark/data_preparation/01-benchmark_pdbs
receptordir=/home/gainza/lpdi_fs/masif_design_paper/repository/masif_seed/masif_seed_search/data/masif_targets_seed_benchmark/data_preparation/01-benchmark_pdbs/
while read ppi_pair;
do
	echo "Starting $ppi_pair"
        outdir=results/$ppi_pair
        PPI_PAIR_ID=$ppi_pair
        PDBID=$(echo $PPI_PAIR_ID | cut -d"_" -f1) 
        CHAIN1=$(echo $PPI_PAIR_ID | cut -d"_" -f2) 
        CHAIN2=$(echo $PPI_PAIR_ID | cut -d"_" -f3) 
        cp $receptordir/$PDBID\_$CHAIN1.pdb $outdir/$PDBID\_$CHAIN1\_m.pdb
        cp $helicesdir/$PDBID\_$CHAIN2.pdb $outdir/$PDBID\_$CHAIN2\_m.pdb
        cd $outdir
        # For the true pair, generate all complexes. 
        ./create.pl zdock_$PDBID\_$CHAIN1\_$PDBID\_$CHAIN2.data $NUM_DECOYS > /dev/null 2>&1
        python3 eval_zrank.py $PPI_PAIR_ID $PDBID\_$CHAIN1\_m.pdb $PDBID\_$CHAIN2\_m.pdb $NUM_DECOYS > results.txt
        break
done < ../benchmark_lists/merged_receptor_helix_list.txt
