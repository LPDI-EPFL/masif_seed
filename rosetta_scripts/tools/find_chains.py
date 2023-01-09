###
# Take a PDB as an input and return a string of the different chain names found in this PDB 
###

import os, sys
import Bio
from Bio.PDB import PDBParser

def find_chains(file):
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', file)
    ch_list=[]
    for model in structure:
        for chain in model:
             ch_list.append(chain.get_id())
    return ch_list

infile=sys.argv[1]
chains=find_chains(infile)
output=chains[0] + ' ' + chains[1]
print(output)

