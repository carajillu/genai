import mdtraj
import scipy
import argparse
import numpy as np
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import sys
def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", type=str, nargs="?", help="System topology (gro or pdb)",required=True)
    parser.add_argument("--traj", type=str, nargs="+", help="Trajectory files (xtc or dcd)",required=True)
    parser.add_argument("--ligand", type=str, nargs="?", help="Ligand name", default=None)
    parser.add_argument("--cutoff", type=float, nargs="?", help="Ligand-Protein interaction cutoff (nm)",default=0.5)
    parser.add_argument("--bsite", type=str, nargs="?", help="Binding site file in gro or pdb format", default=None)
    args = parser.parse_args()
    return args

def concat_traj(topology_path, traj_paths):
    """
    Load multiple dcd/xtc files with the same topology
    and concatenate them into a single trajectory
    """
    traj = mdtraj.load(traj_paths, top=topology_path)
    traj=traj.superpose(reference=traj[0],atom_indices=traj.topology.select("backbone"))
    return traj

def get_ligand_indices(traj, ligand_name):
    """
    Get the indices of the ligand atoms
    """
    ligand_indices = traj.topology.select(f"resname {ligand_name}")
    ligand_indices = [int(i) for i in ligand_indices if traj.topology.atom(i).element != mdtraj.element.hydrogen]
    return ligand_indices

def get_bsite_fromlig(traj, ligand_indices, cutoff):
    """
    For each frame in traj, get the indices of the heavy atoms that 
    are within cutoff nm of any atoms in ligand_indices, but are not
    in ligand_indices themselves. Return a list with unique indices
    """
    ligand_name=traj.topology.atom(ligand_indices[0]).residue.name
    frame_neighbors = mdtraj.compute_neighbors(traj, cutoff, ligand_indices)
    all_neighbors = [int(i) for i in np.concatenate(frame_neighbors, axis=0) if traj.topology.atom(i).residue.name!=ligand_name]
    all_neighbors = list(set(all_neighbors))
    return all_neighbors

def save_bsite_gro(traj, binding_site_indices):
    """
    get a trajectory of the ligand in the binding site.
    Output a gro file with the first frame and a 
    xtc file with the whole trajectory
    """
    slice=binding_site_indices
    newtraj=traj.atom_slice(slice)
    newtraj[0].save_gro("bsite.gro")
    return
    
if __name__ == "__main__":
    args = parse()
    print(args)
    traj = concat_traj(args.top, args.traj)
    ligand_indices = get_ligand_indices(traj, args.ligand)
    binding_site_indices = get_bsite_fromlig(traj, ligand_indices, args.cutoff)   
    save_bsite_gro(traj, binding_site_indices)