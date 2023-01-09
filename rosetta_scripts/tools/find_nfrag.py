import sys

def find_nfrag(seed_pdb):
    input_pdb = open(seed_pdb.strip())
    i=0
    for line in input_pdb.readlines():
        if(line[0:3]=='TER'):
            i+=1
    return i

array=sys.argv[1]
seed_file='./cropseed_'+array+'.pdb'
nfrag=find_nfrag(seed_file)
print(nfrag)
