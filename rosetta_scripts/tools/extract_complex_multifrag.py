###
# Take a seed:target complex PDB as an input and an array number, and will output separated complex in two different PDBs with a name based on the array number.
# Seed PDBs are cropped and composed only of residues found in beta sheet regions (based on DSSP)
###

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import argparse
import os
import pandas as pd
import numpy as np
import re
import sys

def line_read_dssp(ding):
    """
    Takes a pdb filepath and returns name and DSSP secondary structure output
    """
    # Break up lines from file to get pdb and path
    pdb_filename = ding.split('/')[-1]
    pdb_filename = pdb_filename.replace("\n", "")
    pdb_id = pdb_filename[0:4]
    pdb_dir = ding.split('/')[:-1]
    pdb_dir = '/'.join([str(n) for n in pdb_dir])

    # Go to directory and load pdb to get DSSP
    os.chdir(pdb_dir)
    p = PDBParser()
    structure = p.get_structure(pdb_id, pdb_filename)
    model = structure[0]
    dssp = DSSP(model, pdb_filename, dssp="/work/upcorreia/bin/sequence/dssp")
    full_dssp = list(dssp)
    sec_str_result = [i[2] for i in full_dssp]
    sec_str_result = ''.join([str(n) for n in sec_str_result])

    return sec_str_result

def find_length( file, chain ):

    list_resi=[]
    input_file = open(file)
    for line in input_file.readlines():
        split_line=line.split()

        if len(split_line) > 3 and split_line[0] == "ATOM" and split_line[4]==chain:

            list_resi.append(split_line[5])
    length = int(list_resi[-1]) - int(list_resi[0]) + 1
    return length

def find_start( file, chain ):

    list_resi=[]
    input_file = open(file)
    for line in input_file.readlines():
        split_line=line.split()

        if len(split_line) > 3 and split_line[0] == "ATOM" and split_line[4]==chain:

            list_resi.append(split_line[5])

    return int(list_resi[0])

base_dir = os.getcwd()
complex_file=sys.argv[1]
array=sys.argv[2]
ch_seed='B'
ch_target='A'

dssp=line_read_dssp(complex_file)
os.chdir(base_dir)
start_seed=find_start(complex_file, ch_seed)
len_seed=find_start(complex_file, ch_seed)
dssp_seed=dssp[start_seed-1:start_seed+len_seed]
print(dssp_seed)

filename_context='context_'+str(array)+'.pdb'
with open(filename_context,'w') as file_context:
    with open(complex_file,'r') as in_pdb:
        for line in in_pdb:
            if (len(line.split())>3):
                if(((line[0:4] == "ATOM") or (line[0:6] == "HETATM")) and (line[21:22] == ch_target)):
                    file_context.write(line)
        file_context.write('TER')

filename_seed='seed_'+str(array)+'.pdb'
previous_dssp=''
with open(filename_seed,'w') as file_seed:
    with open(complex_file,'r') as in_pdb:
        for line in in_pdb:
            if((line[0:4] == "ATOM") and (line[21:22] == ch_seed)):
                resid=int(line[22:26].strip(' '))
                if(dssp[resid-1]=='E'):
                    file_seed.write(line)
                if(previous_dssp=='E' and dssp[resid-1]!='E'):
                    file_seed.write('TER\n')
                previous_dssp=dssp[resid-1]
