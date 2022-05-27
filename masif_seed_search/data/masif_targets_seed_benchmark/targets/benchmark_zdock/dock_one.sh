PDBID_CHAIN1_CHAIN2=$1
PDBID=$(echo $PDBID_CHAIN1_CHAIN2 | cut -d"_" -f1)
CHAIN1=$(echo $PDBID_CHAIN1_CHAIN2 | cut -d"_" -f2)
CHAIN2=$(echo $PDBID_CHAIN1_CHAIN2 | cut -d"_" -f3)

outdir=03-results/$PDBID_CHAIN1_CHAIN2/
mkdir -p $outdir
# Copy target pdb, with surface marked and blocked residues.
cp 02-zdock_marked_blocked_pdbs/$PDBID\_$CHAIN1\_m_bl.pdb $outdir

# Also copy target pdb, with surface marked but without blocked residues.
cp 01-zdock_marked/$PDBID\_$CHAIN1\_m.pdb $outdir

# Copy ligand pdb of the co-crystal binder, with surface marked
cp 01-zdock_marked/$PDBID\_$CHAIN2\_m.pdb $outdir

cp submit_dock_one.slurm $outdir

# Copy zdock files 
echo $PDBID\_$CHAIN2 > $outdir/input_list.txt
cat ../benchmark_lists/peptide_decoys_random_1000.txt >> $outdir/input_list.txt
cd $outdir 
ln -s ../../zdock
ln -s ../../uniCHARMM
sbatch submit_dock_one.slurm $PDBID $CHAIN1

