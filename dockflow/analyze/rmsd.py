import mdtraj as md
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from meeko import MoleculePreparation, PDBQTMolecule
import numpy as np
import os
def atoms_to_rdkit_mol(atoms_array) -> Chem.Mol:
    mol = Chem.RWMol()
    conformer = Chem.Conformer(len(atoms_array))

    for i, atom_data in enumerate(atoms_array):
        symbol = atom_data['name'].strip()
        try:
            atom = Chem.Atom(symbol)
        except:
            raise ValueError(f"Invalid atom symbol: {symbol}")

        mol_idx = mol.AddAtom(atom)

        # Set 3D coordinates
        x, y, z = atom_data['xyz']
        conformer.SetAtomPosition(mol_idx, Point3D(float(x), float(y), float(z)))

        # Set optional partial charge as a property
        mol.GetAtomWithIdx(mol_idx).SetDoubleProp("PartialCharge", float(atom_data['partial_charges']))

    mol.AddConformer(conformer, assignId=True)
    return mol.GetMol()

def reorder_atoms_to_match_reference(ref_traj: md.Trajectory, target_traj: md.Trajectory) -> md.Trajectory:
    """
    Reorders atoms in the target PDB to match the atom order of the reference PDB.
    
    Parameters
    ----------
    ref_pdb : str
        Path to the reference PDB file (with desired atom order).
    target_pdb : str
        Path to the PDB file to be reordered (same atoms, different order).
    
    Returns
    -------
    md.Trajectory
        MDTraj trajectory object of the target, reordered to match the reference.
    """
    # Load both trajectories
    new_top=md.Topology()
    chain=new_top.add_chain("A")
    residue=new_top.add_residue(ref_traj.topology.residue(0).name,chain)
    residue.resSeq=ref_traj.topology.residue(0).resSeq
    for atom_ref in ref_traj.topology.atoms:
        for atom_target in target_traj.topology.atoms:
            if atom_ref.name == atom_target.name:
                new_atom=new_top.add_atom(atom_ref.name,atom_ref.element,residue)
                new_atom.serial=atom_ref.serial
                break
    new_traj=md.Trajectory(xyz=target_traj.xyz,topology=new_top)
    return new_traj



def compute_rmsds_against_ref(pdb_path, pdbqt_path,n_poses=1000):
    # Load reference ligand (PDB) with MDTraj and remove Hs
    ref_traj = md.load(pdb_path)
    ref_traj = ref_traj.atom_slice([a.index for a in ref_traj.topology.atoms if a.element.symbol != 'H'])
    ref_traj.save_pdb("ref.pdb")

    # Load PDBQT file and extract all poses (MODEL blocks)
    poses = PDBQTMolecule.from_file(pdbqt_path)
    rmsd_list = []
    for pose in poses:
        mol = atoms_to_rdkit_mol(pose.atoms())
        # Remove Hs from RDKit mol
        mol = Chem.RemoveAllHs(mol)
        #Chem.SanitizeMol(mol)
        Chem.MolToPDBFile(mol, "mol.pdb")
        print(os.getcwd())
        new_traj = reorder_atoms_to_match_reference(ref_traj, md.load("mol.pdb"))
        new_traj.save_pdb("new_traj.pdb")
        # Compute RMSD
        rmsd = md.rmsd(new_traj, ref_traj)[0]
        rmsd_list.append(rmsd)

    return rmsd_list
