import open3d as o3d
import copy 
from IPython.core.debugger import set_trace
import numpy as np
from simple_mesh import Simple_mesh
from pathlib import Path
import scipy.sparse as spio
import pymesh
from Bio.PDB import PDBParser, PDBIO, Selection
import tensorflow as tf
from tensorflow import keras
import os


def load_protein_pcd(full_pdb_id, chain_number, paths, flipped_features=False, read_mesh=True):
    """
    Returns the protein point cloud, list of descriptors, interface predictions
    and patch coordinates given full pdb name (str), chain number (int) and if descriptors
    should be flipped or not (bool).
    """

    full_pdb_id_split = full_pdb_id.split('_')
    pdb_id = full_pdb_id_split[0]
    chain_ids = full_pdb_id_split[1:]
    pdb_chain = chain_ids[chain_number - 1]

    pdb_pcd = o3d.read_point_cloud(
        str(Path(paths['surf_dir']) /
            '{}_{}.ply'.format(
            pdb_id,
            pdb_chain)))

    if read_mesh == True:
        pdb_mesh = o3d.read_triangle_mesh(
            str(Path(paths['surf_dir']) /
            '{}_{}.ply'.format(
            pdb_id,
            pdb_chain))) 

    pdb_iface = np.squeeze(np.load(
        Path(paths['iface_dir']) /
        'pred_{}_{}.npy'.format(
            pdb_id,
            pdb_chain)))

    if flipped_features:
        pdb_desc_sc05 = np.load(
            Path(paths['desc_dir']) /
            full_pdb_id /
            'p{}_desc_flipped.npy'.format(chain_number))
    else:
        pdb_desc_sc05 = np.load(
            Path(paths['desc_dir']) /
            full_pdb_id /
            'p{}_desc_straight.npy'.format(chain_number))


    pdb_desc = np.array([pdb_desc_sc05,])

    if read_mesh:
        return pdb_pcd, pdb_desc, pdb_iface, pdb_mesh
    else:
        return pdb_pcd, pdb_desc, pdb_iface 

def get_patch_geo(
        pcd,
        patch_coords,
        center,
        descriptors,
        outward_shift=0.25,
        flip_normals=False):
    """
    Returns a patch from a point cloud pcd with center point center (int),
    based on geodesic distances from patch coords and corresponding Feature descriptors.
    """
    patch_idxs = patch_coords[center]
    patch_pts = np.asarray(pcd.points)[patch_idxs, :]
    patch_nrmls = np.asarray(pcd.normals)[patch_idxs, :]
    patch_pts = patch_pts + outward_shift * patch_nrmls
    if flip_normals:
        patch_nrmls = -patch_nrmls

    patch = o3d.PointCloud()
    patch.points = o3d.Vector3dVector(patch_pts)
    patch.normals = o3d.Vector3dVector(patch_nrmls)
    patch_descs = [o3d.Feature(), o3d.Feature(), o3d.Feature()]
    patch_descs[0].data = descriptors[0,patch_idxs, :].T
    return patch, patch_descs, patch_idxs


def get_patch_coords(top_dir, pdb, pid, cv=None):
    """ 
    Load precomputed patch coordinates.
    """
    if cv is None:
        patch_coords = np.load(os.path.join(
            top_dir, pdb, pid+'_list_indices.npy'), allow_pickle=True)
        cv = np.arange(0, patch_coords.shape[0])
    else:
        patch_coords = np.load(os.path.join(
            top_dir, pdb, pid+'_list_indices.npy'), allow_pickle=True)[cv]
    patch_coords = {key: patch_coords[ii] for ii, key in enumerate(cv)}
    return patch_coords 

# Get the vertex indices of target sites.
def get_target_vix(pc, iface, num_sites=1, selected_vertices=None):
    target_vertices = []
    iface = np.copy(iface)
    if selected_vertices is None:
        selected_vertices = np.arange(len(iface))
    for site_ix in range(num_sites):
        iface_patch_vals = []
        # Go through each patch.
        best_ix = -1
        best_val = float("-inf")
        best_neigh = []
        for ii in selected_vertices:
            neigh = pc[ii]
            val = np.mean(iface[neigh])
            if val > best_val:
                best_ix = ii
                best_val = val
                best_neigh = neigh

        # Now that a site has been identified, clear the iface values so that the same site is not picked again.
        iface[best_neigh] = 0
        target_vertices.append(best_ix)

    return target_vertices

def test_alignments(transformation, source_structure, target_pcd_tree, interface_dist=10.0):
    """
        test_alignments - evaluate quality of the alignments.
    """
    structure_coords = np.array([atom.get_coord() for atom in source_structure.get_atoms()])
    structure_coord_pcd = o3d.PointCloud()
    structure_coord_pcd.points = o3d.Vector3dVector(structure_coords)
    structure_coord_pcd.transform(transformation)
        
    d_nn_interface, i_nn_interface = target_pcd_tree.query(np.asarray(structure_coord_pcd.points),k=1,distance_upper_bound=interface_dist)
    interface_atoms = np.where(d_nn_interface<=interface_dist)[0]
    rmsd = np.sqrt(np.mean(np.square(np.linalg.norm(structure_coords[interface_atoms,:]-np.asarray(structure_coord_pcd.points)[interface_atoms,:],axis=1))))
    return rmsd

def match_descriptors(directory_list, pids, target_desc, params):
    """ 
        Match descriptors to the target descriptor.
    """

    all_matched_names = []
    all_matched_vix = []
    all_matched_desc_dist = []
    count_proteins = 0
    count_descriptors = 0

    for ppi_pair_id in directory_list:
        if 'seed_pdb_list' in params and ppi_pair_id not in params['seed_pdb_list']:
            continue

        if '.npy' in ppi_pair_id or '.txt' in ppi_pair_id:
            continue

        mydescdir = os.path.join(params['seed_desc_dir'], ppi_pair_id)
        for pid in pids:
            try:
                fields = ppi_pair_id.split('_')
                if pid == 'p1':
                    pdb_chain_id = fields[0]+'_'+fields[1]
                elif pid == 'p2':
                    pdb_chain_id = fields[0]+'_'+fields[2]
                iface = np.load(params['seed_iface_dir']+'/pred_'+pdb_chain_id+'.npy')[0]
                descs = np.load(mydescdir+'/'+pid+'_desc_straight.npy')
            except:
                continue
            #print(pdb_chain_id)
            name = (ppi_pair_id, pid)
            count_proteins += 1
            count_descriptors += len(iface)

            diff = np.sqrt(np.sum(np.square(descs - target_desc), axis=1))

            true_iface = np.where(iface > params['iface_cutoff'])[0]
            near_points = np.where(diff < params['desc_dist_cutoff'])[0]

            selected = np.intersect1d(true_iface, near_points)
            if len(selected > 0):
                all_matched_names.append([name]*len(selected))
                all_matched_vix.append(selected)
                all_matched_desc_dist.append(diff[selected])
                #print('Matched {}'.format(ppi_pair_id))
                #print('Scores: {} {}'.format(iface[selected], diff[selected]))
            if count_proteins % 1000 == 0:
                print('Compared target with {} \'fragments\' from {} proteins'.format(count_descriptors, count_proteins))


    # Flatten the matches and convert to a dictionary 
    if len(all_matched_names) == 0: 
        print('Iterated over {} fragments from {} proteins; matched 0 based on descriptor similarity.'.format(count_descriptors, count_proteins))
        print('No descriptors were matched. Check your parameters.')
        return {}
    matched_names = np.concatenate(all_matched_names, axis=0)
    matched_vix = np.concatenate(all_matched_vix, axis=0)
    print('Iterated over {} fragments from {} proteins; matched {} based on descriptor similarity.'.format(count_descriptors, count_proteins, len(matched_vix)))
    matched_desc_dist = np.concatenate(all_matched_desc_dist, axis=0)
    matched_dict = {}
    for name_ix, name in enumerate(matched_names):
        name = (name[0], name[1])
        if name not in matched_dict:
            matched_dict[name] = []
        matched_dict[name].append(matched_vix[name_ix])
            
    return matched_dict

def count_clashes(transformation, source_surface_vertices, source_structure, \
        target_ca_pcd_tree, target_pcd_tree, radius=2.0):

    structure_ca_coords = np.array([atom.get_coord() for atom in source_structure.get_atoms() if atom.get_id() == 'CA'])
    structure_ca_coord_pcd = o3d.PointCloud()
    structure_ca_coord_pcd.points = o3d.Vector3dVector(structure_ca_coords)
    structure_ca_coord_pcd_notTransformed = copy.deepcopy(structure_ca_coord_pcd)
    structure_ca_coord_pcd.transform(transformation)
        
    structure_atoms = [atom for atom in source_structure.get_atoms() if not atom.get_name().startswith('H')]
    structure_coords = np.array([atom.get_coord() for atom in structure_atoms])
    structure_coord_pcd = o3d.PointCloud()
    structure_coord_pcd.points = o3d.Vector3dVector(structure_coords)
    structure_coord_pcd_notTransformed = copy.deepcopy(structure_coord_pcd)
    structure_coord_pcd.transform(transformation)

    # Compute the distance between every source_surface_vertices and every target vertex.
    d_vi_at, i_vi_at = target_pcd_tree.query(source_surface_vertices, k=1)

    # Transform the input structure. - This is a bit of a bad programming practice.
    for ix, v in enumerate(structure_coord_pcd.points):
        structure_atoms[ix].set_coord(v)

    # Invoke cKDTree with (structure_coord_pcd)

    d_nn_ca, i_nn_ca = target_pcd_tree.query(np.asarray(structure_ca_coord_pcd.points),k=1,distance_upper_bound=radius)
    d_nn, i_nn = target_pcd_tree.query(np.asarray(structure_coord_pcd.points),k=1,distance_upper_bound=radius)
    clashing_ca = np.sum(d_nn_ca<=radius)
    clashing = np.sum(d_nn<=radius)
    
    return clashing_ca,clashing,d_vi_at


def multidock(source_pcd, source_patch_coords, source_descs, 
            cand_pts, target_pcd, target_descs,
               params):
    ransac_radius=params['ransac_radius'] 
    ransac_iter=params['ransac_iter']
    all_results = []
    all_source_patch = []
    all_source_scores = []
    all_source_desc = []
    all_source_idx = []
    for pt in cand_pts:
        source_patch, source_patch_descs, source_patch_idx = get_patch_geo(
            source_pcd, source_patch_coords, pt, source_descs, outward_shift=params['surface_outward_shift'])
           
        result = o3d.registration_ransac_based_on_feature_matching(
            source_patch, target_pcd, source_patch_descs[0], target_descs[0],
            ransac_radius,
            o3d.TransformationEstimationPointToPoint(False), 3,
            [o3d.CorrespondenceCheckerBasedOnEdgeLength(0.9),
             o3d.CorrespondenceCheckerBasedOnDistance(1.0),
             o3d.CorrespondenceCheckerBasedOnNormal(np.pi/2)],
            o3d.RANSACConvergenceCriteria(ransac_iter, 500)
        )
        result_icp = o3d.registration_icp(source_patch, target_pcd,
                    1.0, result.transformation, o3d.TransformationEstimationPointToPlane())

        source_patch.transform(result_icp.transformation)
        all_results.append(result_icp)
        all_source_patch.append(source_patch)
        all_source_desc.append(source_patch_descs)
        all_source_idx.append(source_patch_idx)

    return all_results, all_source_patch, all_source_desc, all_source_idx 


def align_and_save(out_filename_base, patch, source_structure, 
                   point_importance, 
                   ):
    
    io = PDBIO()
    # Remove hydrogens.
    for atom in Selection.unfold_entities(source_structure, 'A'):
        if atom.get_name().startswith('H'):
            parent = atom.get_parent()
            parent.detach_child(atom.get_id())
    io.set_structure(source_structure)
    io.save(out_filename_base+'.pdb')
    # Save patch
    mesh = Simple_mesh(np.asarray(patch.points))
    mesh.set_attribute('vertex_charge', point_importance)
    mesh.save_mesh(out_filename_base+'_patch.ply')
    

# Compute different types of scores: 
# -- Compute the neural network score
def compute_nn_score(target_pcd, source_pcd, corr, 
        target_desc, source_desc, 
        target_ckdtree, nn_score, 
        vertex_to_atom_dist, nn_score_cutoff=1.0):


    # Compute nn scores 
    # Compute all points correspondences and distances for nn
    d, r = target_ckdtree.query(source_pcd.points)

    # r: for every point in source, what is its correspondence in target
    # feat0 - distance between poitnts
    distance = d # distances between points 
    distance[distance<0.5] = 0.5
    distance = 1.0/distance

    # feat1: descriptors
    desc_dist = target_desc[0].data[:,r] - source_desc[0].data[:,:]
    desc_dist = np.sqrt(np.sum(np.square(desc_dist.T), axis=1))
    desc_dist = 1.0/desc_dist
    # desc dist score: sum over those within 1.5A
    neigh = np.where(d < 1.5)[0]
    desc_dist_score = np.sum(np.square(desc_dist[neigh]))

    # feat3: normal dot product
    n1 = np.asarray(source_pcd.normals)
    n2 = np.asarray(target_pcd.normals)[r]
    normal_dp = np.multiply(n1, n2).sum(1)

    # feat9 - vertex to atom distance
    vert_atom_dist = vertex_to_atom_dist
    vert_atom_dist[vert_atom_dist < 0.25] = 0.25
    vert_atom_dist = 1.0/vert_atom_dist

    features = np.vstack([distance, desc_dist, normal_dp, vert_atom_dist]).T

    nn_score_pred, point_importance = \
                    nn_score.eval_model(
#                            features, nn_score_cutoff
                            features, 0.9 
                            )

    ret = (np.array([nn_score_pred[0][0], desc_dist_score]).T, point_importance)
    return ret 

def align_protein(name, \
        target_patch, \
        target_patch_descs, \
        target_ckdtree, \
        target_ca_pcd_tree, \
        target_pcd_tree, \
        source_paths, \
        matched_dict, \
        nn_score, \
        site_outdir, \
        params):
    ppi_pair_id = name[0]
    pid = name[1]
    pdb = ppi_pair_id.split('_')[0]

    # PDB Parser
    parser = PDBParser()

    if pid == 'p1':
        chain = ppi_pair_id.split('_')[1]
        chain_number = 1
    else: 
        chain = ppi_pair_id.split('_')[2]
        chain_number = 2
        
    # Load source ply file, coords, and descriptors.
    
    #print('{}'.format(pdb+'_'+chain))
    source_pcd, source_desc, source_iface = load_protein_pcd(ppi_pair_id, chain_number, source_paths, flipped_features=False, read_mesh=False)

    # Get coordinates for all matched vertices.
    source_vix = matched_dict[name]
    source_coord =  get_patch_coords(params['seed_precomp_dir'], ppi_pair_id, pid, cv=source_vix)
    
    # Perform all alignments to target. 
    all_results, all_source_patch, all_source_patch_desc, all_source_idx = multidock(
            source_pcd, source_coord, source_desc,
            source_vix, target_patch, target_patch_descs, 
            params
            ) 

    # Score the results using a 'lightweight' scoring function.
    all_source_scores = [None]*len(source_vix)
    for viii in range(len(source_vix)):
        if len(all_results[viii].correspondence_set)/float(len(np.asarray(all_source_patch[viii].points))) < 0.3:
            # Ignore those with poor fitness.
            all_source_scores[viii] = ([0,0], np.zeros_like(all_source_idx[viii]))
        else:
            # Compute the distance between every source_surface_vertices and every target vertex.
            d_vi_at, _= target_pcd_tree.query(np.asarray(all_source_patch[viii].points), k=1)

            all_source_scores[viii] = compute_nn_score(target_patch, all_source_patch[viii], np.asarray(all_results[viii].correspondence_set),
                 target_patch_descs, all_source_patch_desc[viii], \
                target_ckdtree, nn_score, d_vi_at, 1.0) # Ignore point importance for speed.
            if len(np.asarray(all_results[viii].correspondence_set)) <= 2.0:
                if all_source_scores[viii][0][0] > 0.9: 
                    set_trace()
#                    set_trace()

    # All_point_importance: impact of each vertex on the score. 
    all_point_importance = [x[1] for x in all_source_scores]
    all_source_scores = [x[0] for x in all_source_scores]
    scores = np.asarray(all_source_scores)
    
    # Filter anything below neural network score cutoff
    top_scorers = np.where(scores[:,0] > params['nn_score_cutoff'])[0]
    
    if len(top_scorers) > 0:
        source_outdir = os.path.join(site_outdir, '{}'.format(ppi_pair_id))
        if not os.path.exists(source_outdir):
            os.makedirs(source_outdir)

        # Load source structure 
        # Perform the transformation on the atoms
        for j in top_scorers:
            res = all_results[j]

            # Load the source pdb structure 
            source_struct = parser.get_structure('{}_{}'.format(pdb,chain), os.path.join(params['seed_pdb_dir'],'{}_{}.pdb'.format(pdb,chain)))
            # Count clashes and align source_structure_cp
            source_structure_cp = copy.deepcopy(source_struct)
            clashing_ca, clashing_total, _= count_clashes(res.transformation, np.asarray(all_source_patch[j].points), source_structure_cp, \
                target_ca_pcd_tree, target_pcd_tree)

            # Check if the number of clashes exceeds the number allowed. 
            if clashing_ca <= params['allowed_CA_clashes'] and clashing_total <= params['allowed_heavy_atom_clashes']:

                # Align and save the pdb + patch 
                out_fn = source_outdir+'/{}_{}_{}'.format(pdb, chain, j)
                align_and_save(out_fn, all_source_patch[j], source_structure_cp, \
                                            point_importance=all_point_importance[j])

                # Recompute the score with clashes.
                print('Selected fragment: {} fragment_id: {} score: {:.4f} desc_dist_score: {:.4f} clashes(CA): {} clashes(total):{}\n'.format(j , ppi_pair_id, scores[j][0], scores[j][1], clashing_ca, clashing_total))

                # Align and save the ply file for convenience.     
                mesh = Simple_mesh()
                mesh.load_mesh(os.path.join(params['seed_surf_dir'],'{}.ply'.format(pdb+'_'+chain)))
                
                # Output the score for convenience. 
                out_score = open(out_fn+'.score', 'w+')
                out_score.write('name: {}, point id: {}, score: {:.4f}, clashing_ca: {}, clashing_heavy: {}, desc_dist_score: {}\n'.format(ppi_pair_id, j, scores[j][0], clashing_ca,clashing_total, scores[j][1]))
                out_score.close()

                source_pcd_copy = copy.deepcopy(source_pcd)
                source_pcd_copy.transform(res.transformation)
                out_vertices = np.asarray(source_pcd_copy.points)
                out_normals = np.asarray(source_pcd_copy.normals)
                mesh.set_attribute('vertex_x', out_vertices[:,0])
                mesh.set_attribute('vertex_y', out_vertices[:,1])
                mesh.set_attribute('vertex_z', out_vertices[:,2])
                mesh.set_attribute('vertex_nx', out_normals[:,0])
                mesh.set_attribute('vertex_ny', out_normals[:,1])
                mesh.set_attribute('vertex_nz', out_normals[:,2])
                mesh.vertices = out_vertices
                mesh.set_attribute('vertex_iface', source_iface) 
                mesh.save_mesh(out_fn+'.ply')
                #save_ply(out_fn+'.ply', out_vertices, mesh.faces, out_normals, charges=mesh.get_attribute('vertex_charge'))
