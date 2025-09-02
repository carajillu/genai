import mdtraj
import numpy as np
import glob
import pandas as pd


def pdbqt_to_mdtraj_trajectory(pdbqt_path: str = "vina_out.pdbqt"):
    """
    Load a multi-pose Vina PDBQT file and return an MDTraj Trajectory with one frame per pose.

    Parameters:
    - pdbqt_path (str): Path to the multi-pose PDBQT file.

    Returns:
    - mdtraj.Trajectory: MDTraj trajectory object (coordinates in nanometers).
    """
    poses = []
    atom_names = []
    elements = []

    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()

    current_coords = []
    element_switch=True
    for line in lines:
        if line.startswith("MODEL"):
            current_coords = []
        elif line.startswith("ENDMDL"):
            element_switch=False
            poses.append(np.array(current_coords, dtype=np.float32))
        elif line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            current_coords.append([x, y, z])
            if element_switch:  # record atom info only once
                atom_name = line[12:16].strip()
                atom_names.append(atom_name)
                elements.append(atom_name)
    #print(elements)
    if not poses:
        raise ValueError("No poses found in PDBQT file.")
    
    # Convert coordinates to nanometers
    xyz = np.array(poses) / 10.0  # shape: (n_frames, n_atoms, 3)

    # Create MDTraj Topology
    topology = mdtraj.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("LIG", chain)
    for name, element in zip(atom_names, elements):
        topology.add_atom(name, mdtraj.element.get_by_symbol(element), residue)

    return mdtraj.Trajectory(xyz=xyz, topology=topology)

def merge_mdtraj_trajectories(traj1, traj2):
    """
    Merge two MDTraj trajectories frame-by-frame into one trajectory
    that contains all atoms from both.

    Parameters:
    - traj1: mdtraj.Trajectory
    - traj2: mdtraj.Trajectory

    Returns:
    - mdtraj.Trajectory: merged trajectory
    """
    if len(traj1) != len(traj2):
        raise ValueError("Trajectories must have the same number of frames.")

    # Combine coordinates
    combined_xyz = np.concatenate([traj1.xyz, traj2.xyz], axis=1)  # shape: (n_frames, n_atoms1 + n_atoms2, 3)

    # Combine topologies
    combined_topology = mdtraj.Topology()
    chain1 = combined_topology.add_chain()
    for res in traj1.topology.residues:
        new_res = combined_topology.add_residue(res.name, chain1)
        for atom in res.atoms:
            combined_topology.add_atom(atom.name, atom.element, new_res)

    chain2 = combined_topology.add_chain()
    for res in traj2.topology.residues:
        new_res = combined_topology.add_residue(res.name, chain2)
        for atom in res.atoms:
            combined_topology.add_atom(atom.name, atom.element, new_res)

    # Create new trajectory
    merged_traj = mdtraj.Trajectory(xyz=combined_xyz, topology=combined_topology)

    return merged_traj

def compute_contact_surface_area_single_pose(
    bsite_traj: mdtraj.Trajectory,
    ligand_traj: mdtraj.Trajectory,
    probe_radius: float = 0.14,
    n_sphere_points: int = 960
) -> float:
    """
    Compute the contact surface area (CSA) between a ligand and selected receptor atoms in a single pose.

    Parameters
    ----------
    bsite_traj : mdtraj.Trajectory
        MDTraj trajectory for the binding site (1 frame).
    ligand_traj : mdtraj.Trajectory
        MDTraj trajectory for the ligand (1 frame).
    probe_radius : float, optional
        Probe radius in nm for SASA calculation (default is 0.14).
    n_sphere_points : int, optional
        Number of points used in Shrake-Rupley algorithm (default is 960).

    Returns
    -------
    float
        Contact surface area in nm^2.
    """
    assert bsite_traj.n_frames == 1 and ligand_traj.n_frames == 1, "Both inputs must be single-frame trajectories."

    # Concatenate ligand and selected receptor atoms
    traj_combined = merge_mdtraj_trajectories(ligand_traj,bsite_traj)

    # Compute SASAs
    sasa_ligand = mdtraj.shrake_rupley(ligand_traj, probe_radius=probe_radius,
                                   n_sphere_points=n_sphere_points, mode='atom').sum()
    sasa_receptor = mdtraj.shrake_rupley(bsite_traj, probe_radius=probe_radius,
                                     n_sphere_points=n_sphere_points, mode='atom').sum()
    sasa_combined = mdtraj.shrake_rupley(traj_combined, probe_radius=probe_radius,
                                     n_sphere_points=n_sphere_points, mode='atom').sum()

    csa = sasa_ligand + sasa_receptor - sasa_combined
    return float(csa)

def compute_vina_out_csa(receptor: str="receptor.pdb", receptor_name: str=None, receptor_selection: str=None, ligand_list: list[str]=["ligand_out.pdbqt"], pose: int=0):
    if receptor_name is None:
       receptor_name = receptor.split("/")[-1].split(".")[0]
    receptor_trj=mdtraj.load(receptor)
    csa_list=[]
    for ligand_str in ligand_list:
        try:
           ligand_trj=pdbqt_to_mdtraj_trajectory(ligand_str)[pose]
           csa_list.append(compute_contact_surface_area_single_pose(receptor_trj,ligand_trj))
        except Exception as e:
            print(f"WARNING for ligand {ligand_str}: {e}")
            csa_list.append(0)
    z=pd.DataFrame()
    z["ligand"]=ligand_list
    z[f"CSA_{receptor_name}"]=csa_list
    return z