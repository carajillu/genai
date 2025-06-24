from pathlib import Path
import numpy as np
import mdtraj as md

def load_pdbqt_as_mdtraj(pdbqt_path: Path) -> md.Trajectory:
    """
    Load a multi-pose PDBQT file and return an MDTraj trajectory.

    Parameters
    ----------
    pdbqt_path : Path
        Path to the .pdbqt file.

    Returns
    -------
    md.Trajectory
        Trajectory with one frame per pose.
    """
    coords = []
    current_pose = []
    atom_names = []
    elements = []

    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                current_pose = []

            elif line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                current_pose.append([x, y, z])

                if not coords:  # only for first model
                    name = line[12:16].strip()
                    atom_names.append(name)
                    try:
                        elements.append(md.core.element.get_by_symbol(name[:2].capitalize()))
                    except KeyError:
                        raise ValueError(f"Unknown element: {name}")

            elif line.startswith("ENDMDL"):
                if current_pose:
                    coords.append(current_pose)
                    current_pose = []

    if current_pose:  # in case last pose lacks ENDMDL
        coords.append(current_pose)

    coords = np.array(coords, dtype=np.float32) / 10.0  # convert from Å to nm

    # Build MDTraj topology
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue("LIG", chain)

    for name, element in zip(atom_names, elements):
        top.add_atom(name, element, residue)

    traj = md.Trajectory(xyz=coords, topology=top)
    return traj
