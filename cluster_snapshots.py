import mdtraj
import scipy
import argparse
import numpy as np
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import sys
import subprocess
def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", type=str, nargs="?", help="System topology (gro or pdb)",required=True)
    parser.add_argument("--traj", type=str, nargs="+", help="Trajectory files (xtc or dcd)",required=True)
    parser.add_argument("--bsite", type=str, nargs="?", help="Binding site file in gro or pdb format", default=None)
    parser.add_argument("--gen-vina-conf", action="store_true", help="Generate a vina conf file", default=False)
    args = parser.parse_args()
    assert args.bsite is not None or args.ligand is not None, "Either ligand or binding site must be provided"
    return args

def get_bsite_fromfile(bsite_path, traj):
    """
    Get the indices of the binding site from a file
    """
    bsite_obj = mdtraj.load(bsite_path)
    bsite_indices=[]
    #This loop is because indices in apo and holo are not the same
    for atom1 in bsite_obj.topology.atoms:
        if atom1.element==mdtraj.element.hydrogen:
            continue
        for atom2 in traj.topology.atoms:
            if atom2.element==mdtraj.element.hydrogen:
                continue
            if atom1.residue.name==atom2.residue.name and atom1.residue.resSeq==atom2.residue.resSeq and atom1.name==atom2.name:
                bsite_indices.append(atom2.index)
    return bsite_indices

def cluster_traj(traj, binding_site_indices):
    """
    Cluster the binding site trajectory (no hydrogens)
    """
    bsite_indices=[ i for i in binding_site_indices if traj.topology.atom(i).element!=mdtraj.element.hydrogen]
    bsite_trj=traj.atom_slice(bsite_indices)
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

    return medoid_indices

def generate_bsite_pdbqt(traj, bsite_indices, medoid_indices):
    """
    Generate a pdbqt file of the binding site (with hydrogens) for each medoid frame
    """
    #bsite_trj=traj.atom_slice(bsite_indices)
    #chop end of chains (vina doesn't understand them)
    uknwnrs=["NTH","CAR","NHI","CGL"]
    sel=[]
    for atom in traj.topology.atoms:
        print(atom.residue.name)
        if atom.residue.name[0:3] not in uknwnrs:
            sel.append(atom.index)
    traj=traj.atom_slice(sel)
    for index in medoid_indices:
        traj[index].save_pdb(f"receptor_{index}.pdb")
        subprocess.run(["mk_prepare_receptor.py", "-i", f"receptor_{index}.pdb", "--write_pdbqt", f"receptor_{index}.pdbqt"])

def generate_vina_conf(traj, bsite_indices, medoid_indices):
    """
    for each medoid frame, find the maximum and minimum coordinates of the binding site
    and use them to set the center and size of the box. The center is the average of the max
    and min coordinates, multiplied by 10 because vina uses Angstroms and the coordinates are in nm.
    The size extends from the centre and is in grid points, which are 1 angstrom apart.
    for each medoid frame, generate a vina conf file called vina_bsite{medoid_index}.conf 
    with the receptor name as bsite{medoid_index}.pdbqt, ligand name as ligand.pdbqt,
    and the center and size of the box.
    The format is as follows:

    receptor = bsite{medoid_index}.pdbqt
    ligand = ligand.pdbqt
    center_x = x
    center_y = y
    center_z = z
    size_x = x
    size_y = y
    size_z = z
    num_modes = 10
    """
    for index in medoid_indices:
        bsite_trj = traj.atom_slice(bsite_indices)
        max_coords = np.array([np.max(traj.xyz[:, :, 0]), np.max(traj.xyz[:, :, 1]), np.max(traj.xyz[:, :, 2])])*10
        min_coords = np.array([np.min(traj.xyz[:, :, 0]), np.min(traj.xyz[:, :, 1]), np.min(traj.xyz[:, :, 2])])*10
        center = (max_coords + min_coords) / 2
        size = (max_coords - center) * 1 # vina uses grid points, which are 1 angstrom apart
        with open(f"vina_bsite_{index}.conf", "w") as f:
            f.write(f"receptor=bsite_{index}.pdbqt\n")
            f.write(f"ligand=ligand.pdbqt\n")
            f.write(f"center_x={center[0]}\ncenter_y={center[1]}\ncenter_z={center[2]}\n")
            f.write(f"size_x={int(size[0])}\nsize_y={int(size[1])}\nsize_z={int(size[2])}\n")
            f.write(f"num_modes=100\n")

if __name__ == "__main__":
    args = parse()
    print(args)
    traj = mdtraj.load(args.traj, top=args.top)
    bsite_indices = get_bsite_fromfile(args.bsite, traj)
    medoid_indices = cluster_traj(traj, bsite_indices)
    if args.gen_vina_conf:
        generate_bsite_pdbqt(traj, bsite_indices, medoid_indices)
        generate_vina_conf(traj, bsite_indices, medoid_indices)