import os, sys, glob, shutil, pickle
import numpy as np
import pandas as pd
import open3d as o3d
import plotly.express as px
import plotly.graph_objects as go

sys.path.append("/scratch/shxiao/masif_seed/masif_seed_search/source")
sys.path.append("/scratch/shxiao/masif_seed/masif/source")
os.environ["masif_target_root"] = "/scratch/shxiao/masif_seed/masif/"
os.environ["masif_seed_search_root"] = "/scratch/shxiao/masif_seed/masif_seed_search/"

from default_config.masif_opts import masif_opts
from matplotlib import pyplot as plt
from scipy.spatial import cKDTree
from Bio.PDB import *
from Bio import SeqUtils
from Bio.PDB import PDBIO, PDBParser, Selection
from Bio.PDB.Polypeptide import is_aa
from tqdm.notebook import tqdm
from scipy.spatial.distance import cdist

sys.path.append('/work/lpdi/users/shxiao/notebooks')
from utils import *



target_pcd, target_desc, target_iface = load_protein_pcd(target_name, target_chain_ix, target_paths, flipped_features=True, read_mesh=False)
source_pcd, source_desc, source_iface = load_protein_pcd(source_name, source_chain_ix, source_paths, flipped_features=False, read_mesh=False)

target_patch, target_patch_descs, target_patch_idx = get_patch_geo(
    target_pcd, target_patch_coords, pt, target_descs, outward_shift=params['surface_outward_shift'], flip_normals=True)

source_patch, source_patch_descs, source_patch_idx = get_patch_geo(
    source_pcd, source_patch_coords, pt, source_descs, outward_shift=params['surface_outward_shift'], flip_normals=False)

result = registration_ransac_based_on_feature_matching(
    source_patch, target_patch, source_patch_descs[0], target_patch_descs[0], True,
    ransac_radius,
    TransformationEstimationPointToPoint(False), 3,
    [CorrespondenceCheckerBasedOnEdgeLength(0.9),
    CorrespondenceCheckerBasedOnDistance(1.0),
    CorrespondenceCheckerBasedOnNormal(np.pi/2)],
    RANSACConvergenceCriteria(max_iteration=params['ransac_iter'], confidence=0.9999),
    seed=42
)

ransac_transformation = result.transformation 

# TODO: there is a potential bug here in benchmark cases only. If a random rotation is not previously applied, 
# there exists a possibility that masif-search/masif-site will get the right patches, but ransac will fail to 
# get a transformation. In these rare cases, icp will start from the ground truth. To get around this, make
# sure a rotation is applied beforehand (for benchmarks at least)
result_icp = registration_icp(source_patch, target_patch,
            1.0, result.transformation, TransformationEstimationPointToPlane())

source_patch.transform(result_icp.transformation)
all_results.append(result_icp)
all_source_patch.append(source_patch)
all_source_desc.append(source_patch_descs)
all_source_idx.append(source_patch_idx)