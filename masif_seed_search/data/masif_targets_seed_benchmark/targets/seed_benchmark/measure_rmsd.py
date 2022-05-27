from Bio.PDB import * 
import numpy as np
from scipy.spatial import cKDTree
import sys 
# Masif-seed helices benchmark
# Measure the RMSD of the ground truth with that of the predicted model in the helices benchmark.

def compute_irmsd(
    receptor_ca_atoms, 
    gt_helix_ca_atoms,
    predicted_helix_ca_atoms,
    interface_dist=10.0, # The interface cutoff to define the interface.
):
    """
    Verify the alignment against the ground truth. 
    """

    
    receptor_ckd = cKDTree(receptor_ca_atoms) 
    dist, res = receptor_ckd.query(gt_helix_ca_atoms)
    iface_ix = np.where(dist < interface_dist)[0]

    rmsd = np.sqrt(
        np.mean(
            np.square(
                np.linalg.norm(
                    gt_helix_ca_atoms[iface_ix]
                    - predicted_helix_ca_atoms[iface_ix],
                    axis=1,
                )
            )
        )
    )
    return rmsd

pdbid_receptor = sys.argv[1]
pdbid_gt_helix = sys.argv[2]
pdbid_predicted_helix = sys.argv[3]

gtdir_helix_tmpl = "/home/gainza/lpdi_fs/masif_design_paper/masif/data/masif_helices/data_preparation/01-benchmark_pdbs/{}.pdb"
gtdir_receptor_tmpl = "{}.pdb"
align_helix_tmpl0 = "site_0/{}/{}.pdb"
align_helix_tmpl1 = "site_1/{}/{}.pdb"

parser= PDBParser()
rec_fn = gtdir_receptor_tmpl.format(pdbid_receptor, pdbid_receptor)
rec_struct = parser.get_structure(rec_fn, rec_fn)
rec_res = Selection.unfold_entities(rec_struct, 'R')
rec_ca_atom_coord = np.asarray([x['CA'].get_coord() for x in rec_res if 'CA' in x])

gt_helix_fn = gtdir_helix_tmpl.format(pdbid_gt_helix)
gt_helix_struct = parser.get_structure(gt_helix_fn, gt_helix_fn)
gt_helix_res = Selection.unfold_entities(gt_helix_struct, 'R')
gt_helix_ca_atom_coord = np.asarray([x['CA'].get_coord() for x in gt_helix_res if 'CA' in x])

predicted_helix_fn0 = align_helix_tmpl0.format(pdbid_gt_helix, pdbid_predicted_helix)
predicted_helix_fn1 = align_helix_tmpl1.format(pdbid_gt_helix, pdbid_predicted_helix)
try:
    predicted_helix_struct = parser.get_structure(predicted_helix_fn0, predicted_helix_fn0)
except:
    predicted_helix_struct = parser.get_structure(predicted_helix_fn1, predicted_helix_fn1)
predicted_helix_res = Selection.unfold_entities(predicted_helix_struct, 'R')
predicted_helix_ca_atom_coord = np.asarray([x['CA'].get_coord() for x in predicted_helix_res if 'CA' in x])

print(compute_irmsd(rec_ca_atom_coord, gt_helix_ca_atom_coord, predicted_helix_ca_atom_coord))
