import os 
import numpy as np
from Bio.PDB import *
from scipy.spatial import cKDTree
import sys
# Compute the results of the masif seed benchmark and the rmsd to the ground truth.
# Measure the RMSD of the ground truth with that of the predicted model in the helices benchmark.
def measure_rmsd(ppi_pair_id, aligned_dirname, aligned_pdb_fn, site):
    gtdir_helix_tmpl = "/home/gainza/lpdi_fs/masif_design_paper/masif/data/masif_helices/data_preparation/01-benchmark_pdbs/{}.pdb"
    gtdir_receptor_tmpl = "{}.pdb"
    predicted_helix_fn = f"{site}/{aligned_dirname}/{aligned_pdb_fn.replace('.score', '.pdb')}"
    pdbid_receptor  = '_'.join([ppi_pair_id.split('_')[0], ppi_pair_id.split('_')[1]])
    pdbid_gt_helix = '_'.join([ppi_pair_id.split('_')[0], ppi_pair_id.split('_')[2]])

    parser= PDBParser()
    rec_fn = gtdir_receptor_tmpl.format(pdbid_receptor)
    rec_struct = parser.get_structure(rec_fn, rec_fn)
    rec_res = Selection.unfold_entities(rec_struct, 'R')
    rec_ca_atom_coord = np.asarray([x['CA'].get_coord() for x in rec_res if 'CA' in x])

    gt_helix_fn = gtdir_helix_tmpl.format(pdbid_gt_helix)
    gt_helix_struct = parser.get_structure(gt_helix_fn, gt_helix_fn)
    gt_helix_res = Selection.unfold_entities(gt_helix_struct, 'R')
    gt_helix_ca_atom_coord = np.asarray([x['CA'].get_coord() for x in gt_helix_res if 'CA' in x])

    predicted_helix_struct = parser.get_structure(predicted_helix_fn, predicted_helix_fn)
    predicted_helix_res = Selection.unfold_entities(predicted_helix_struct, 'R')
    predicted_helix_ca_atom_coord = np.asarray([x['CA'].get_coord() for x in predicted_helix_res if 'CA' in x])

    return compute_irmsd(rec_ca_atom_coord, gt_helix_ca_atom_coord, predicted_helix_ca_atom_coord)

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

ppi_pair_id = sys.argv[1]
pdbid = ppi_pair_id.split('_')[0]
all_scores = []
for site in ['site_0', 'site_1']:
    if os.path.exists(site):
        for dirname in os.listdir(site):
            curpath = os.path.join(site, dirname)
            if not dirname.endswith('.vert'):
                for fn in os.listdir(curpath): 
                    if fn.endswith('.score'): 
                        fullfn = os.path.join(curpath, fn)
                        with open(fullfn) as f: 
                            result = f.readlines()[0]
                        score = float(result.split(',')[2].split()[1])
                        if pdbid in dirname: 
                            irmsd = measure_rmsd(ppi_pair_id, dirname, fn, site)
                            all_scores.append((dirname, score, irmsd))
                        else: 
                            all_scores.append((dirname, score, 1e8))

all_scores = sorted(all_scores, key=lambda x: x[1], reverse=True)
for i in range(len(all_scores)): 
    if pdbid in all_scores[i][0] and all_scores[i][2] < 3.0:
        print(f'{pdbid} rank: {i}, score: {all_scores[i][1]}, rmsd: {all_scores[i][2]}')
        sys.exit(0)
print(f'{pdbid} N/A')

