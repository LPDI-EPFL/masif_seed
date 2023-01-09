# Rosetta Scripts documentation
Two groups of Rosetta scripts have been used in this publication : one for refining the seed provided by MaSIF_seed and one for grafting the seeds on protein scaffolds.
## Seed refinement
### Definition
This scripts aims to:
1) Take as an input the seeds provided by MaSIF
2) Put the seed PDB in complex with the target and relax the pose with a FastRelax protocol
3) Re-design the seed interface with a FastDesign protocol and a penalty for buried unsatisfied polar atoms in the Rosetta scoring function
### Use
1) Place the seed PDBs provided by MaSIF and the target PDB in the `input` directory
2) Create a file `list.txt` that contains two columns, for the seed filename and the target filename respectively:
```
./input/seedPDB1.pdb ./input/targetPDB.pdb
./input/seedPDB2.pdb ./input/targetPDB.pdb
...
```
3) Modify the submission file `submitter.slurm` according to your cluster settings, and notably:
   * The path to the Rosetta binary files in `ROSETTA_BIN`
   * The number of arrays (= number of seeds) in `#SBATCH --array=...`
4) Submit the job to your cluster using `submitter.slurm`
### Output
The script will output the relaxed seed-target complex in the `relax` folder and the redesigned seed in complex with the target in the `output` folder. Some common metrics such as shape complementarity, number of hydrogen bonds, number of buried unsatisfied polar atoms and computed binding energy could be found in the score file `score.sc`. We recommend to select the best 30-50 seeds prior grafting. 
### Reference
   * Coventry, B. & Baker, D. Protein sequence optimization with a pairwise decomposable penalty for buried unsatisfied hydrogen bonds. PLOS Comput. Biol. 17, e1008061 (2021).

## Seed grafting
### Definition 
This scripts aims to:
1) Take as an input a selection of refined seeds
2) Rename the seed with a unique ID number and...
   * Alpha helices seeds : Define the hotspot residues for each seed
   * Beta sheet seeds: Crop the loop region, define the hotspot residues for each seed and crop the side-chain without any contact
3) Graft each seed in parallel on a suitable scaffold originating from a scaffold database (Double job paralellization) using the MotifGraft mover
4) Refine the scaffold with two rouds of design and minimization with a penalty for buried unsatisfied polar atoms in the Rosetta scoring function
### Use
1) Place the selected seed PDB (recommended : 30-50 seeds) in the `input folder`
2) Modify the submission file `run_master_motif_graft.slurm` according to your cluster settings, and notably :
   * The number of arrays (= number of seeds) in `#SBATCH --array=...`
3) Similarly, modify the submission file `run_motif_graft_varhotspots.slurm` according to your cluster settings, and notably :
   * The path `in_dir` to the database of scaffolds used for grafting
   * The number of arrays `#SBATCH --array=...` for paralellization in regards to the number of splitfiles that composed the scaffold database
   * The path to the Rosetta binary files in `ROSETTA_BIN`
4) Run the grafting protocol using `run_protocol.sh`
### Output
The script will output several designs per seed that sucessfully passed the grafting criteria. Results folder will be split per seed ID number and per array number. The python script `gather_score.py` in the `tool` folder can help to gather all separate score files in one unique document. Some common metrics such as shape complementarity, number of hydrogen bonds, number of buried unsatisfied polar atoms and computed binding energy can be used for selecting the final designs. 
### References
   * Silva, D., Correia, B.E., and Procko, E. (2016) Motif-driven Design of Protein-Protein Interactions. Methods Mol. Biol. 1414:285-304
