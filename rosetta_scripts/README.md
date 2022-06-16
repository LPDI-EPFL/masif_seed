## Rosetta Scripts documentation
Two groups of Rosetta scripts have been used in this publication : one for refining the seed provided by MaSIF_seed and one for grafting the seed on a protein scaffold.
### Seed refinement
This scripts aims to :
1) Take as an input the seeds provided by MaSIF
2) Put the seed PDB in complex with the target and relax the pose 
3) Re-design with interface with a penalty for buried unsatisfied polar atoms in the Rosetta scoring function
