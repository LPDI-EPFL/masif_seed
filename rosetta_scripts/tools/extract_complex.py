###
# Take a seed:target complex PDB as an input and an array number, and will output separated complex in two different PDBs with a name based on the array number. 
###

import os, sys
import Bio
from Bio.PDB import PDBParser

#base_dir = os.getcwd()
#os.chdir(base_dir)

def find_chains(file):
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', file)
    ch_list=[]
    for model in structure:
        for chain in model:
             ch_list.append(chain.get_id())
    return ch_list

complex_file=sys.argv[1]
array=sys.argv[2]
#complex_file='1DCI002_A_23_3Q0Y_B.pdb'
#array=1

### If/else order to be changed if the seed is the first or second chain 
i=0
for chain in find_chains(complex_file):
    if (i==0):
        filename_context='context_'+str(array)+'.pdb'
        with open(filename_context,'w') as file_context:
            with open(complex_file,'r') as in_pdb:
                for line in in_pdb:
                    if (len(line.split())>3):
                        if(((line[0:4] == "ATOM") or (line[0:6] == "HETATM")) and (line[21:22] == chain)):
                            file_context.write(line)
            file_context.write('TER')
        i=1
    else:
        filename_seed='seed_'+str(array)+'.pdb'
        with open(filename_seed,'w') as file_seed:
            with open(complex_file,'r') as in_pdb:
                for line in in_pdb:
                    if (len(line.split())>3):
                        if((line[0:4] == "ATOM") and (line[21:22] == chain)):
                            file_seed.write(line)
            file_seed.write('TER')
        
