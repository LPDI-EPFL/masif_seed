###
# Take a seed:target complex PDB and both chain names as an input, and return a string containing all the hotspots making contact with the target (<2.5A).
# Formating of the string is compatible with the 'hotspots' argument of the MotifGraft mover of RosettaScript.
###

import os,sys,string,gzip
import math

#base_dir = os.getcwd()
#os.chdir(base_dir)

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

def find_min_dist(input_pdb, resid, chain_target):

    dist=1000
    res_id=0
    
    input_pdb1 = open(input_pdb.strip())
    
    
    for line in input_pdb1.readlines():

        if (len(line.split())>3):
            if((line[0:4] == "ATOM") and (line[21:22] != chain_target) and (line[22:26].strip(' ') == str(resid))):
                cs=[float(line[30:38].strip(' ')), float(line[38:46].strip(' ')), float(line[46:54].strip(' '))]
                input_pdb2 = open(input_pdb.strip())
                for line2 in input_pdb2.readlines():
                    if (len(line2.split())>3):
                        if(((line2[0:4] == "ATOM") or line2[0:6] == "HETATM")and (line2[21:22] == chain_target)):
                            ct=[float(line2[30:38].strip(' ')), float(line2[38:46].strip(' ')), float(line2[46:54].strip(' '))]
                            new_dist=math.sqrt(pow((cs[0]-ct[0]),2)+pow((cs[1]-ct[1]),2)+pow((cs[2]-ct[2]),2))
                        
                            if new_dist<dist:
                                dist=new_dist
    return dist

infile = sys.argv[1]
ch_seed = sys.argv[2]
ch_target = sys.argv[3]

len_seed=find_length(infile,ch_seed)
start=find_start(infile,ch_seed)
outstr=''

for i in range(start,start+len_seed,1):
    dist=find_min_dist(infile,i,ch_target)
    if(dist<2.5):
        outstr+=str(i-start+1)
        outstr+=':'
outstr=outstr.strip(':')
print(outstr)
