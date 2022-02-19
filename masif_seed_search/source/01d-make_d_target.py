#!/usr/bin/python
import numpy as np
import os
import Bio
import shutil
from Bio.PDB import * 
import sys
import importlib
from IPython.core.debugger import set_trace

# Local includes
from default_config.masif_opts import masif_opts
from triangulation.computeMSMS import computeMSMS
from triangulation.fixmesh import fix_mesh
import pymesh
from input_output.extractHelix import extractHelix
from input_output.save_ply import save_ply
from input_output.read_ply import read_ply
from input_output.protonate import protonate
from triangulation.computeHydrophobicity import computeHydrophobicity
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from triangulation.computeAPBS import computeAPBS
from triangulation.compute_normal import compute_normal

if len(sys.argv) <= 1: 
    print("Usage: {config} "+sys.argv[0]+" PDBID_A")
    print("A or AB are the chains to include in this surface.")
    sys.exit(1)

pdbid_chain=sys.argv[1]

# First flip the surface file.
# Open the L-surface file
mymesh = pymesh.load_mesh(os.path.join(masif_opts['ply_chain_dir'], pdbid_chain+'.ply'))
# Flip the x coordinate. 
verts = np.copy(mymesh.vertices)
verts[:,0] = -verts[:,0]
# It is fundamental to also flip the x coordinate of the normals.
vertex_nx = - mymesh.get_attribute('vertex_nx')

out_mesh = pymesh.form_mesh(verts, mymesh.faces)
out_mesh.add_attribute('vertex_nx')
out_mesh.set_attribute('vertex_nx', vertex_nx)
out_mesh.add_attribute('vertex_ny')
out_mesh.set_attribute('vertex_ny', mymesh.get_attribute('vertex_ny'))
out_mesh.add_attribute('vertex_nz')
out_mesh.set_attribute('vertex_nz', mymesh.get_attribute('vertex_nz'))

out_mesh.add_attribute('vertex_iface')
out_mesh.set_attribute('vertex_iface', mymesh.get_attribute('vertex_iface'))
out_mesh.add_attribute('vertex_hbond')
out_mesh.set_attribute('vertex_hbond', mymesh.get_attribute('vertex_hbond'))
out_mesh.add_attribute('vertex_hphob')
out_mesh.set_attribute('vertex_hphob', mymesh.get_attribute('vertex_hphob'))
out_mesh.add_attribute('vertex_charge')
out_mesh.set_attribute('vertex_charge', mymesh.get_attribute('vertex_charge'))

outfilename = os.path.join(masif_opts['ply_chain_dir'], 'd'+pdbid_chain+'.ply')
pymesh.save_mesh(
    outfilename, out_mesh, *out_mesh.get_attribute_names(), use_float=True, ascii=True
)

# Now flip the atom coordinates.
parser = PDBParser()
in_pdbfn = os.path.join(masif_opts['pdb_chain_dir'], pdbid_chain+'.pdb')
struct = parser.get_structure(in_pdbfn, in_pdbfn)

atoms = Selection.unfold_entities(struct, 'A')
for atom in atoms: 
    coord = atom.get_coord()
    coord[0] = -coord[0]
    atom.set_coord(coord)

out_pdbfn = os.path.join(masif_opts['pdb_chain_dir'], 'd'+pdbid_chain+'.pdb')
pdbio = PDBIO()
pdbio.set_structure(struct)
pdbio.save(out_pdbfn)
#masif_opts['pdb_chain_dir']

