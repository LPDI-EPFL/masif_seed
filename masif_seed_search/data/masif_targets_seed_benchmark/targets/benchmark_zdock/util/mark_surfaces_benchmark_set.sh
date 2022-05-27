# Generate the pdbs in the format that zdock needs them.
masif_root=$(git rev-parse --show-toplevel)
# The route to your zdock directory
#zdock_dir=$masif_root/ext_programs/zdock3.0.2_mac_intel/
zdock_dir=/home/gainza/lpdi_fs/bin/ZDock/zdock3.0.2_linux_x64

# The benchmark receptor list
receptor_benchmark_list=/home/gainza/lpdi_fs/masif_design_paper/masif/masif_seed_search/data/masif_targets/targets/benchmark_lists/merged_receptor_helix_list.txt

# Random peptides list
random_peptides_list=/home/gainza/lpdi_fs/masif_design_paper/masif/masif_seed_search/data/masif_targets/targets/benchmark_lists/peptide_decoys_random_1000.txt

# Location of receptor pdbs. 
RECEPTORPDBLOC=/home/gainza/lpdi_fs/masif_design_paper/masif/masif_seed_search/data/masif_targets/data_preparation/01-benchmark_pdbs

HELIXPDBLOC=/home/gainza/lpdi_fs/masif_design_paper/masif/data/masif_peptides/data_preparation/01-benchmark_pdbs

mkdir -p output_tmp/

OUTDIR='01-zdock_marked/'
mkdir -p $OUTDIR

while read PDBID_CHAIN1_CHAIN2;
do
  PDBID=$(echo $PDBID_CHAIN1_CHAIN2| cut -d"_" -f1)
  CHAIN1=$(echo $PDBID_CHAIN1_CHAIN2| cut -d"_" -f2)
  CHAIN2=$(echo $PDBID_CHAIN1_CHAIN2| cut -d"_" -f3)

  PDB1=$PDBID\_$CHAIN1\.pdb
  # Remove all hydrogens, as reduce hydrogens are not recognized by marksur. 
  reduce -Trim $RECEPTORPDBLOC/$PDB1 > output_tmp/$PDB1
  # Mark the surface residues for zdock, save in temporary dir
  $zdock_dir/mark_sur output_tmp/$PDB1 $OUTDIR/$PDBID\_$CHAIN1\_m.pdb

  PDB2=$PDBID\_$CHAIN2\.pdb
#  # Remove all hydrogens, as reduce hydrogens are not recognized by marksur. 
  reduce -Trim $HELIXPDBLOC/$PDB2 > output_tmp/$PDB2
#  # Mark the surface residues for zdock, save in temporary dir
  $zdock_dir/mark_sur output_tmp/$PDB2 $OUTDIR/$PDBID\_$CHAIN2\_m.pdb

done < $receptor_benchmark_list

while read PDBID_CHAIN1;
do
  PDBID=$(echo $PDBID_CHAIN1 | cut -d"_" -f1)
  CHAIN1=$(echo $PDBID_CHAIN1 | cut -d"_" -f2)

  PDB1=$PDBID\_$CHAIN1\.pdb
  # Remove all hydrogens, as reduce hydrogens are not recognized by marksur. 
  reduce -Trim $HELIXPDBLOC/$PDB1 > output_tmp/$PDB1
  # Mark the surface residues for zdock, save in temporary dir
  $zdock_dir/mark_sur output_tmp/$PDB1 $OUTDIR/$PDBID\_$CHAIN1\_m.pdb

done < $random_peptides_list
