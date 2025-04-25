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
    assert args.bsite is not None or args.ligand is not None, "Either ligand or binding site must be provided"
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
    frame_neighbors = mdtraj.compute_neighbors(traj, cutoff, ligand_indices)
    all_neighbors = [int(i) for i in np.concatenate(frame_neighbors, axis=0) if (i not in ligand_indices) and (traj.topology.atom(i).element != mdtraj.element.hydrogen)]
    all_neighbors = list(set(all_neighbors))
    return all_neighbors

def get_bsite_fromfile(bsite_path, traj):
    """
    Get the indices of the binding site from a file
    """
    bsite_obj = mdtraj.load(bsite_path)
    bsite_indices=[]
    #This loop is because indices in apo and holo are not the same
    for atom1 in bsite_obj.topology.atoms:
        for atom2 in traj.topology.atoms:
            #if atom1==atom2:
            if atom1.residue.name==atom2.residue.name and atom1.residue.resSeq==atom2.residue.resSeq and atom1.name==atom2.name:
                bsite_indices.append(atom2.index)
    return bsite_indices

def get_bsite_gro(traj, binding_site_indices):
    """
    get a trajectory of the ligand in the binding site.
    Output a gro file with the first frame and a 
    xtc file with the whole trajectory
    """
    slice=binding_site_indices
    newtraj=traj.atom_slice(slice)
    newtraj[0].save_gro("bsite.gro")
    return

def cluster_traj(traj, binding_site_indices):
    """
    Cluster the binding site trajectory
    """
    bsite_trj=traj.atom_slice(binding_site_indices)
    n_frames, n_atoms, _ = bsite_trj.xyz.shape
    X = bsite_trj.xyz.reshape(n_frames, n_atoms*3)
    #run K-Means
    k = 5   # number of clusters you want
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(X)
    labels = kmeans.labels_
    print(np.unique(labels))

    #find a “medoid” (the most central frame) for each cluster
    medoid_indices = {}
    for cid in np.unique(labels):
        if cid == -1:
            continue
        members = np.where(labels == cid)[0]
        # compute pairwise RMSD within this cluster
        sub = bsite_trj[members]
        D = np.empty((len(members), len(members)))
        for i, idx in enumerate(members):
            D[i] = mdtraj.rmsd(sub, sub, i)
        # pick the member whose sum of distances is minimal
        medoid = members[np.argmin(D.sum(axis=1))]
        medoid_indices[cid] = medoid
    print(f"K-Means found {len(medoid_indices)} clusters")
    print("Medoid frame of each cluster:", medoid_indices)

    #calculate the RMSD matrix of all the medoid frames
    rmsd_matrix = np.empty((len(medoid_indices), len(medoid_indices)))
    for i, idx in enumerate(medoid_indices):
        rmsd_matrix[i] = mdtraj.rmsd(bsite_trj[idx], bsite_trj[medoid_indices[i]])
    print(f"Max RMSD: {np.max(rmsd_matrix)}")
    print(f"Min RMSD: {np.min(rmsd_matrix)}")

    for index in medoid_indices:
        bsite_trj[index].save_gro(f"bsite_{index}.gro")
    #bsite_trj[np.array(list(medoid_indices.values()))].save_xtc("medoids.xtc")

    return medoid_indices
if __name__ == "__main__":
    args = parse()
    print(args)
    traj = concat_traj(args.top, args.traj)
    if args.bsite is not None:
        binding_site_indices = get_bsite_fromfile(args.bsite, traj)
        medoid_indices = cluster_traj(traj, binding_site_indices)
        for index in medoid_indices:
            traj[index].save_gro(f"trj_{index}.gro")
            #traj[medoid_indices].save_xtc("medoids.xtc")
    else:
        ligand_indices = get_ligand_indices(traj, args.ligand)
        binding_site_indices = get_bsite_fromlig(traj, ligand_indices, args.cutoff)   
        get_bsite_gro(traj, binding_site_indices)