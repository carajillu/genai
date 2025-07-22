import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import mdtraj as md
import numpy as np
from meeko.preparation import MoleculePreparation
from meeko import PDBQTWriterLegacy
import os
import subprocess

def smiles_to_mdtraj(smiles: str, name: str = "ligand") -> md.Trajectory:
    """
    Converts a SMILES string into an MDTraj Trajectory with hydrogens and 3D coordinates.

    Parameters
    ----------
    smiles : str
        SMILES string for the ligand.
    name : str
        Molecule name (used for atom naming and metadata).

    Returns
    -------
    md.Trajectory
        A single-frame MDTraj trajectory containing the prepared ligand.
    """
    # Parse SMILES and add hydrogens
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if success != 0:
        raise ValueError(f"Embedding failed for ligand: {smiles}")

    # Geometry optimisation (optional but recommended)
    AllChem.UFFOptimizeMolecule(mol)

    # Extract coordinates
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)/10) for i in range(mol.GetNumAtoms())], dtype=np.float32)

    # Build topology from RDKit
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue(name, chain)
    for i, atom in enumerate(mol.GetAtoms()):
        top.add_atom(atom.GetSymbol(), md.element.Element.getBySymbol(atom.GetSymbol()), residue)

    # Create MDTraj trajectory
    traj = md.Trajectory(xyz=coords[None, :, :], topology=top)

    return traj

def smiles_to_pdbqt(smiles: str, output_path: str, name: str = "ligand"):
    """
    Converts a SMILES string into a PDBQT file using RDKit + Meeko (Meeko v0.5+ compatible).

    Parameters
    ----------
    smiles : str
        The SMILES string of the molecule.
    output_path : str
        Path to write the .pdbqt file.
    name : str
        Name of the ligand (used in metadata).
    """
    # Prepare RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
        raise ValueError(f"Could not embed 3D structure for: {smiles}")
    AllChem.UFFOptimizeMoleculeConfs(mol) # seems to optimize better than UFFOptimizeMolecule

    mol.SetProp("_Name", name)

    # Use Meeko to prepare and export PDBQT
    mk_prep = MoleculePreparation()
    molsetup=mk_prep.prepare(mol)[0]
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
    with open(output_path, "w") as f:
        f.write(pdbqt_string)
    return

def sdf_to_pdbqt(sdf_path, pdbqt_path=None):
    """
    Convert a single-molecule SDF file to a PDBQT file using Meeko.

    Parameters:
    - sdf_path (str): Path to the input .sdf file
    - pdbqt_path (str, optional): Output path for the .pdbqt file.
                                  If None, replaces extension with .pdbqt

    Returns:
    - str: Path to the generated .pdbqt file
    """
    # Load molecule with RDKit
    mol = Chem.SDMolSupplier(sdf_path, sanitize=False)[0]
    if mol is None:
        raise ValueError(f"Could not read molecule from {sdf_path}")

    # Prepare with Meeko
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()

    # Determine output path
    if pdbqt_path is None:
        pdbqt_path = os.path.splitext(sdf_path)[0] + ".pdbqt"

    # Write to file
    with open(pdbqt_path, "w") as f:
        f.write(pdbqt_string)

    return pdbqt_path

def prepare_library_from_smiles(smiles_file: str, output_dir: str = "pdbqt_library"):
    """
    Prepare a library of ligands from a SMILES file and save them as PDBQT files.
    The library needs to be formatted as a csv file where the separators are commas 
    and the smiles column is named either 'smiles' or 'SMILES'. 

    Parameters
    ----------
    smiles_file : str
        Path to the input SMILES file.
    output_dir : str
        Directory to save the output PDBQT files.

    Returns
    -------
    str
        Path to the directory containing the PDBQT files.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Read SMILES file
    smiles_df = pd.read_csv(smiles_file)
    assert 'smiles' in smiles_df.columns or 'SMILES' in smiles_df.columns, \
        "SMILES file must contain a column named 'smiles' or 'SMILES'."
    smiles_list = smiles_df['smiles'] if 'smiles' in smiles_df.columns else smiles_df['SMILES']

    for i, smiles in enumerate(smiles_list):
        smiles = smiles.strip()
        if not smiles:
            continue
        output_path = os.path.join(output_dir, f"ligand_{i+1}.pdbqt")
        smiles_to_pdbqt(smiles, output_path)

    return os.path.abspath(output_dir)

def pdbqt_to_sdf(pdbqt_path: str, output_path: str, explicit_H: bool = True):
    """
    Takes as input a multi-frame pdbqt file output by vina and:
    - converts it to a PDB file to add bond orders (with standard obabel command)
    - converts the PDB file to an output sdf file (with explicit hydrogens if specified)

    Returns
    -------
    str
        Path to the generated sdf file.
    """
   
    return output_path