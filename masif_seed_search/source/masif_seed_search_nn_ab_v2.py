#!/usr/bin/env python
#!/usr/bin/env python
import os, sys, dask, time, pymesh, shutil, importlib, argparse
import numpy as np

from scipy.spatial import cKDTree
from IPython.core.debugger import set_trace
from default_config.masif_opts import masif_opts
from simple_mesh import Simple_mesh
from Bio.PDB import PDBParser
from tqdm import tqdm

from open3d import *

from alignment_evaluation_nn import AlignmentEvaluationNN
from alignment_utils_ab_v2 import *

def main():
    parser = argparse.ArgumentParser(description='MaSIF Pipeline for aligning Antibody/Nanobody paratopes to target epitopes.')
    parser.add_argument('custom_params_fn', type=str, help='Custom parameters filename')
    parser.add_argument('target_name', type=str, help='Target protein name')
    parser.add_argument('target_pid', default='p1', type=str, help='Target protein pid; can be either "p1" or "p2"')
    parser.add_argument('--split_seed_list', type=str, help='Split seed matches list')
    parser.add_argument('--site_selection_mode', type=str, help='Can be set to either "dist" or "iface"; if "dist", will select the top N patches by distance to the target; if "iface", will select the top N patches by interface score.')
    
    args = parser.parse_args()

    # Parameters with databases used, etc for seed search. 
    custom_params_fn = args.custom_params_fn
    custom_params_obj = importlib.import_module(custom_params_fn, package=None)
    params = custom_params_obj.params

    # Load target patches.
    target_pid = args.target_pid
    target_full_name = args.target_name
    if target_pid == 'p1':
        target_name = target_full_name.split('_')[0] + '_' + target_full_name.split('_')[1]
        target_chain_ix = 1
    else:
        target_name = target_full_name.split('_')[0] + '_' + target_full_name.split('_')[2]
        target_chain_ix = 2

    if args.split_seed_list:
        # Use a split of the seed matches. 
        seed_list = open(args.split_seed_list)
        seed_ppi_pair_ids  = [x.rstrip() for x in seed_list.readlines()]
    else:
        seed_ppi_pair_ids  = np.array(os.listdir(params['seed_desc_dir']))

    # Initialize two neural networks - one that does not account for atomic clashes (initial filter) and one with clashes. 
    nn_score_atomic = AlignmentEvaluationNN(params['nn_score_atomic_fn'], selected_features=[0,1,2,9], max_npoints=params['max_npoints']) ## Slightly slower but more accurate.
    nn_score_atomic.restore_model()

    # Go through every 12A patch in the target protein -- get a sorted least in order of the highest iface mean in the patch
    target_ply_fn = os.path.join(params['target_ply_iface_dir'], target_name+'.ply')
    mymesh = pymesh.load_mesh(target_ply_fn)
    iface = mymesh.get_attribute('vertex_iface')
    target_coord = get_patch_coords(params['target_precomp_dir'], target_full_name, target_pid)

    # Define target and source paths (for interface scores, descriptors, ply files)
    target_paths = {}
    target_paths['surf_dir'] = params['target_surf_dir'] 
    target_paths['iface_dir'] = params['target_iface_dir'] 
    target_paths['desc_dir'] = params['target_desc_dir'] 

    source_paths = {}
    source_paths['surf_dir'] = params['seed_surf_dir'] 
    source_paths['iface_dir'] = params['seed_iface_dir'] 
    source_paths['desc_dir'] = params['seed_desc_dir'] 

    # Load the target point cloud, descriptors, interface and mesh.
    target_pcd, target_desc, target_iface, target_mesh = load_protein_pcd(target_full_name, target_chain_ix,
                                                                        target_paths, flipped_features=True,
                                                                        read_mesh=True)

    # Open the pdb structure of the target, load into point clouds for fast access.
    parser = PDBParser()
    target_pdb_path = os.path.join(params['target_pdb_dir'],'{}.pdb'.format(target_name))
    target_struct = parser.get_structure(target_pdb_path, target_pdb_path)
    target_atom_coords = [atom.get_coord() for atom in target_struct.get_atoms() if not atom.get_name().startswith('H') ]
    target_ca_coords = [atom.get_coord() for atom in target_struct.get_atoms() if atom.get_id() == 'CA']
    # Create kdtree search trees (for fast comparision).
    target_ca_pcd_tree = cKDTree(np.array(target_ca_coords))
    target_pcd_tree = cKDTree(np.array(target_atom_coords))

    # If a specific residue is selected, then go after that residue 
    if 'target_residue' in params:
        # Use the tuple for biopython: (' ', resid, ' ')
        print(f"Now finding points within {params['target_residue']['cutoff']}A of residue {params['target_residue']['resid']} {params['target_residue']['atom_id']} atom of {target_name.split('_')[0] + '_' + target_name.split('_')[1]}...")
        target_resid = params['target_residue']['resid']
        target_chain = params['target_residue']['chain']
        target_cutoff = params['target_residue']['cutoff']
        target_atom_id = params['target_residue']['atom_id']
        coord = target_struct[0][target_chain][target_resid][target_atom_id].get_coord()
        dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
        if args.site_selection_mode == 'dist':
            target_vertices = dists.argsort()[:params['num_sites']] # Just take the closest point to the residue.
        elif args.site_selection_mode == 'iface':
            neigh_indices = np.where(dists<target_cutoff)[0]
            # Get a target vertex for every target site.
            target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'], selected_vertices=neigh_indices)
            # Add the closest point to the residue as the first point.
            target_vertices = np.insert(target_vertices, 0, dists.argmin(), axis=0)
        else:
            raise ValueError(f"site_selection_mode {args.site_selection_mode} not recognized. Can only be 'dist' or 'iface'.")
        print(f"Point indices: {', '.join([str(t) for t in target_vertices])}")
    elif 'target_sites' in params:
        target_vertices = params['target_sites']
        # coord = mymesh.vertices[params['target_site']]
        # dists = np.sqrt(np.sum(np.square(mymesh.vertices - coord), axis=1))
        # if args.site_selection_mode == 'dist':
        #     target_vertices = dists.argsort()[:params['num_sites']] # Just take the closest point to the residue.
        # elif args.site_selection_mode == 'iface':
        #     neigh_indices = np.where(dists<target_cutoff)[0]
        #     # Get a target vertex for every target site.
        #     target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'], selected_vertices=neigh_indices)
        #     # Add the closest point to the residue as the first point.
        #     target_vertices = np.insert(target_vertices, 0, dists.argmin(), axis=0)
        # else:
        #     raise ValueError(f"site_selection_mode {args.site_selection_mode} not recognized. Can only be 'dist' or 'iface'.")
        print(f"Point indices: {', '.join([str(t) for t in target_vertices])}")
    else:
        # Get a target vertex for every target site.
        target_vertices = get_target_vix(target_coord, iface, num_sites=params['num_sites'])

    outdir = params['out_dir_template'].format(target_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Copy the pdb structure and the ply file of the target 
    shutil.copy(target_pdb_path, outdir)
    shutil.copy(target_ply_fn, outdir)

    # Go through every target site in the target
    if 'selected_site_ixs' in params:
        site_ixs = params['selected_site_ixs']
        site_vixs = [vix for ix,vix in enumerate(target_vertices) if ix in site_ixs]
    else:
        site_ixs = [ix for ix,vix in enumerate(target_vertices)]
        site_vixs = target_vertices

    # Go through every selected site 
    for site_ix, site_vix in zip(site_ixs, site_vixs):
        if 'target_residue' in params:
            residue_name = ''.join([str(x).strip(" ") for x in params['target_residue']['resid']])
            print('Starting site near residue {}\n'.format(residue_name))
            site_outdir = os.path.join(outdir, 'res_{}_site_{}'.format(residue_name, site_vix))
        else:
            print('Starting site {}\n'.format(site_vix))
            site_outdir = os.path.join(outdir, 'res_unk_site_{}'.format(site_vix))
        os.makedirs(site_outdir, exist_ok=True)
        # Get the geodesic patch and descriptor patch for each target patch
        target_patch, target_patch_descs, target_patch_idx = \
                    get_patch_geo(target_pcd, target_coord, site_vix,\
                            target_desc, flip_normals=True, outward_shift=params['surface_outward_shift'])

        # Make a ckdtree with the target vertices.
        target_ckdtree = cKDTree(target_patch.points)

        # Write out the patch itself.  
        out_patch = open(site_outdir+'/target.vert', 'w+')
        for point in target_patch.points: 
            out_patch.write('{}, {}, {}\n'.format(point[0], point[1], point[2]))
        out_patch.close()

        # Match the top descriptors in the database based on descriptor distance.
        print(''.join(['-']*50))
        print('Starting to match target descriptor to descriptors from {} proteins; this may take a while.'.format(len(seed_ppi_pair_ids))) 
        matched_dict = match_descriptors(seed_ppi_pair_ids, ['p1'], target_desc[0][site_vix], params)
    
        if len(matched_dict.keys())==0:
            continue 

        print(''.join(['-']*50))
        print("Second stage of MaSIF seed search: each matched descriptor is aligned and scored; this may take a while..")
        count_matched_fragments = 0 
        for ix, name in enumerate(matched_dict.keys()):
            print(f'\nStarting to align {len(matched_dict[name])} fragments from {name[0]}...')
            align_protein(
                name,
                target_patch, 
                target_patch_descs, 
                target_ckdtree, 
                target_ca_pcd_tree, 
                target_pcd_tree, 
                source_paths, 
                matched_dict,
                nn_score_atomic, 
                site_outdir, 
                params
            )
            if (ix+1) %1000 == 0:
                print('So far, MaSIF ahs aligned {} fragments from {} proteins.'.format(count_matched_fragments,ix+1))
            count_matched_fragments+= len(matched_dict[name]) 

if __name__ == "__main__":
    main()
