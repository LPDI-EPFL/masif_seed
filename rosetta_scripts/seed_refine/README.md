## Documentation on using the seed refinement protocol
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
