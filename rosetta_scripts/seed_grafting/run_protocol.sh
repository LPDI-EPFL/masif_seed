bash ../tools/seednum.sh
ls ./input/*.pdb > list.txt
sbatch run_master_motif_graft.slurm
