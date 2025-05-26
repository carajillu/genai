import argparse
import mdtraj
import numpy as np
import pandas as pd

def parse():
    parser = argparse.ArgumentParser(description="Compute contact surface area between two selections in an MD trajectory.")
    parser.add_argument("--topology", required=True, help="Topology file (e.g., pdb, gro)")
    parser.add_argument("--trajectory", required=False, default=None, help="Trajectory file (e.g., xtc, dcd)")
    parser.add_argument("--selection_1", required=True, help="MDTraj selection string for the first group")
    parser.add_argument("--selection_2", required=True, nargs="+", help="MDTraj selection string for the second group")
    parser.add_argument("--timestep", required=False, default=0.002, type=float, help="Time step in picoseconds (ps)")
    parser.add_argument("--output", default="contact_surface_area.csv", help="Output CSV file to write results")

    return parser.parse_args()

def compute_sasa(traj, atom_indices, probe_radius=0.14, n_sphere_points=960):
    """
    Compute solvent-accessible surface area (SASA) for given atoms in each frame.
    """
    sasa = mdtraj.shrake_rupley(traj.atom_slice(atom_indices), probe_radius=probe_radius, n_sphere_points=n_sphere_points, mode='atom')
    sasa_total=np.sum(sasa, axis=1)  # Sum over atoms per frame
    return sasa_total


def main():
    args = parse()

    print("Loading trajectory...")

    if args.trajectory is None:
        traj = mdtraj.load(args.topology)
    else:
        traj = mdtraj.load(args.trajectory, top=args.topology)

    
    print(f"Selecting atoms for {args.selection_1} ...")
    atoms_1 = traj.topology.select(args.selection_1)
    assert len(atoms_1) > 0, "No atoms found for selection 1"
    print(f"Number of atoms in selection 1: {len(atoms_1)}")
    for atom in atoms_1:
        print(traj.topology.atom(atom),traj.topology.atom(atom).residue.index)
    sasa_1 = compute_sasa(traj, atoms_1)

    contacts=pd.DataFrame()
    contacts["Time_ps"]=np.arange(len(sasa_1)) * args.timestep
    
    for selection in args.selection_2:
        print(f"Selecting atoms for {selection} ...")
        atoms_2 = traj.topology.select(selection)
        assert len(atoms_2) > 0, "No atoms found for selection 2"
        print(f"Number of atoms in selection 2: {len(atoms_2)}")
        for atom in atoms_2:
            print(traj.topology.atom(atom),traj.topology.atom(atom).residue.index)
        sasa_2 = compute_sasa(traj, atoms_2)
        sasa_combined = compute_sasa(traj, np.concatenate([atoms_1, atoms_2]))
        contact_surface = (sasa_1 + sasa_2) - sasa_combined
        contacts[selection.replace(" ","_")]=contact_surface

    print(f"Writing results to {args.output}...")
    contacts.to_csv(args.output, index=False)

    print("Done.")

if __name__ == "__main__":
    main()
