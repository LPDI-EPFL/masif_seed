mkdir -p analysis_fixed_size/
mkdir -p analysis_fixed_size/out_pdb/
mkdir -p analysis_fixed_size/out_data/
find out_peptides/ -name "*.score" -printf "%f\n" | sed 's/.score//g' > analysis_fixed_size/list_out_peptides.txt
num_lines=$(wc -l analysis_fixed_size/list_out_peptides.txt | cut -d" " -f1)
sed -i "s/#SBATCH --array=.*/#SBATCH --array=1-$num_lines/" align_all_to_all_fixed_length.slurm
