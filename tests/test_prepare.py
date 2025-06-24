import pytest
import os
import tempfile
from rdkit import Chem
import mdtraj as md
from dockflow.prepare.ligands import smiles_to_mdtraj, smiles_to_pdbqt  # Replace with the actual module path

def test_smiles_to_mdtraj():
    """Test that a valid SMILES returns a correct MDTraj Trajectory."""
    smiles = "CCO"  # Ethanol
    name = "ethanol"

    traj = smiles_to_mdtraj(smiles, name)
    # Build RDKit molecule to compare atom count
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    n_atoms = mol.GetNumAtoms()

    assert isinstance(traj, md.Trajectory)
    assert traj.n_frames == 1
    assert traj.n_atoms == n_atoms
    assert traj.xyz.shape == (1, n_atoms, 3)



def test_smiles_to_pdbqt():
    """Test that a PDBQT file is created and contains expected content."""
    smiles = "CCO"  # Ethanol

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, "ligand.pdbqt")
        smiles_to_pdbqt(smiles, out_path, name="ethanol")
        # Check that file exists
        assert os.path.exists(out_path)
       # Check that file is not empty
        with open(out_path, "r") as f:
            content = f.read()
        
        assert len(content) > 0
        assert "ROOT" in content or "REMARK" in content or "ATOM" in content
