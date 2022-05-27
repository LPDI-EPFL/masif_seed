import numpy as np
import pymesh
import sys
import time
import shutil
import subprocess
import os
#import pyflann
from IPython.core.debugger import set_trace


# Script to block all atoms that are not in the target site-- used to block everything in ZDOCK outside the target site.
# This is done by adding a 19 to column 55-56 for atoms outside target.
# Input in the form 4QVF_A_B, where 4QVF is the pdbid and A, B are the chains. The first chain ('A') is blocked.
target_pdb_pair = sys.argv[1]
target_pdb = target_pdb_pair.split('_')[0]+'_'+target_pdb_pair.split('_')[1]
outdir = '02-zdock_marked_blocked_pdbs/'
# masif site prediction dir (to find the exact target site).
masif_site_ply_dir = '../../output/all_feat_3l_seed_benchmark/pred_surfaces/'
precomputation_dir = '../../data_preparation/04b-precomputation_12A/precomputation/'

# Read the surface for this target.
# Use Masif-site to identify the target site. 
mesh = pymesh.load_mesh(os.path.join(masif_site_ply_dir, target_pdb+'.ply'))
iface = mesh.get_attribute('vertex_iface')
# Load the patch indices
patch_indices = np.load(os.path.join(precomputation_dir, target_pdb, 'p1_list_indices.npy'), allow_pickle=True)
# Find the center of the best patch. 
best_patch = -1
best_score = -1
best_neighbors = -1
for ix, patch in enumerate(patch_indices): 
    score = np.mean(iface[patch])
    if score > best_score:
        best_score = score
        best_patch = ix
        best_neighbors = patch
        
    # This is the predicted interface center.
    center_vertex = best_patch
    center_point = mesh.vertices[center_vertex]

# Read the neighbors for this center point.
neigh = best_neighbors
print(neigh)
# Read the vertices
v_t = mesh.vertices

# Open the pdb file line by line. 
pdbid = target_pdb_pair.split('_')[0]
t_chain = target_pdb_pair.split('_')[1]
pdb_file_t = os.path.join('01-zdock_marked/','{}_{}_m.pdb'.format(pdbid,t_chain))
all_lines = open(pdb_file_t, 'r').readlines()
# Store line_ix and coord for all lines
interface_line_ix = []
target_atom_coords = []
for ix, line in enumerate(all_lines):
    if line.startswith("ATOM"):
        interface_line_ix.append(ix)
        coordx = float(line[30:38])
        coordy = float(line[38:46])
        coordz = float(line[46:54])

        target_atom_coords.append((coordx, coordy, coordz))

target_atom_coords = np.array(target_atom_coords)

# Find the atom coords closest to each point in the neighborhood of the interface
neigh_atom_ix = set()
for point in neigh:
    dists = np.sqrt(np.sum(np.square(target_atom_coords - v_t[point]),axis=1 ))
    closest_at_iii = np.argmin(dists)
    neigh_atom_ix.add(interface_line_ix[closest_at_iii])

# Block every line that is not in the interface, by adding a 19 to columns 55-56
outfile = open(outdir+'{}_{}_m_bl.pdb'.format(pdbid,t_chain),'w')
for ix, line in enumerate(all_lines):
    if line.startswith("ATOM"):
        # Block every line whose index is not in neigh_atom_ix
        if ix not in neigh_atom_ix: 
            listline = list(line)
            listline[55] = '1'
            listline[56] = '9'
            line = ''.join(listline)
        else:
            print(line)
        outfile.write(line)
        
