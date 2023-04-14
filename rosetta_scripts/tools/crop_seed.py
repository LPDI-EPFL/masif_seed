import os,sys,string,gzip
import math

#base_dir = os.getcwd()
#os.chdir(base_dir)

def find_length( file, chain, frag_id ):

    list_resi=[]
    input_file = open(file)
    current_frag=1
    for line in input_file.readlines():
        split_line=line.split()

        if len(split_line) > 3 and split_line[0] == "ATOM" and split_line[4]==chain and current_frag==frag_id:
            list_resi.append(split_line[5])
        if line[0:3]=='TER':
            current_frag+=1
    length = int(list_resi[-1]) - int(list_resi[0]) + 1
    start = int(list_resi[0])
    return start,length

def find_min_dist(context_pdb, seed_pdb, resid, chain_target):

    dist=1000
    res_id=0
    
    input_pdb1 = open(seed_pdb.strip())
    
    
    for line in input_pdb1.readlines():
        if (len(line.split())>3):
            if((line[0:4] == "ATOM") and (line[21:22] != chain_target) and (line[22:26].strip(' ') == str(resid)) and (line[12:16].strip(' ') not in ['C','N','HA','CA','O','H'])):
                cs=[float(line[30:38].strip(' ')), float(line[38:46].strip(' ')), float(line[46:54].strip(' '))]
                input_pdb2 = open(context_pdb.strip())
                for line2 in input_pdb2.readlines():
                    if (len(line2.split())>3):
                        if(((line2[0:4] == "ATOM") or line2[0:6] == "HETATM")and (line2[21:22] == chain_target)):
                            ct=[float(line2[30:38].strip(' ')), float(line2[38:46].strip(' ')), float(line2[46:54].strip(' '))]
                            new_dist=math.sqrt(pow((cs[0]-ct[0]),2)+pow((cs[1]-ct[1]),2)+pow((cs[2]-ct[2]),2))
                        
                            if new_dist<dist:
                                dist=new_dist
    return dist

def find_nfrag(seed_pdb):
    input_pdb = open(seed_pdb.strip())
    i=0
    for line in input_pdb.readlines():
        if(line[0:3]=='TER'):
            i+=1
    return i

array = sys.argv[1]
ch_seed = 'B'
ch_target = 'A'
seed_file='./seed_'+array+'.pdb'
cropseed_file='./cropseed_'+array+'.pdb'
context_file='./context_'+array+'.pdb'
#infile=r"3DH7002_A_7_3OSK_V.pdb"
#ch_seed='A'
#ch_target='B'

outstr=''

with open(cropseed_file,'w') as outf:
    for j in range (1, find_nfrag(seed_file)+1,1):
        start,len_seed=find_length(seed_file,ch_seed, j)
        hotspots=[]
        for i in range(start,start+len_seed,1):
            dist=find_min_dist(context_file,seed_file,i,ch_target)
            print(i,dist)
            if(dist<2.5):
                hotspots.append(i)
        if(len(hotspots)==1):
            with open(seed_file,'r') as in_pdb:
                for line in in_pdb:
                    if(line[0:4] == "ATOM"):
                        resid=int(line[22:26].strip(' '))
                        if((resid==hotspots[0]-1) or (resid==hotspots[0]) or (resid==hotspots[0]+1)):
                            outf.write(line)
            outf.write('TER\n')                
        if(len(hotspots)>1):
            with open(seed_file,'r') as in_pdb:
                for line in in_pdb:
                    if(line[0:4] == "ATOM"):
                        resid=int(line[22:26].strip(' '))
                        if((resid>=hotspots[0]) and (resid<=hotspots[-1])):
                            outf.write(line)
            outf.write('TER\n')
            
