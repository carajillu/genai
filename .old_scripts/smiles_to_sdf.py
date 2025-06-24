import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import subprocess
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import sys

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_str", type=str, nargs="?", help="Smiles string", default=None)
    parser.add_argument("--smiles_path", type=str, nargs="?", help="Smiles file", default=None)
    args = parser.parse_args()
    assert args.smiles_str is not None or args.smiles_path is not None, "Smiles string or file must be provided"
    return args

def smiles_to_3D(smiles, mol_name):
    print(smiles)
    """
    Convert a smiles string to an sdf file
    """
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp('_Name', mol_name)
    molh=Chem.AddHs(mol)
    AllChem.EmbedMolecule(molh,AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(molh,maxIters=1000)
    return molh

def sdf_to_pdbqt(sdf_path,pdbqt_path):
    for mol in Chem.SDMolSupplier(sdf_path,removeHs=False):
        mk_prep = MoleculePreparation()
        molsetup=mk_prep.prepare(mol)[0]
        pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
        with open(pdbqt_path, "w") as f:
            f.write(pdbqt_string)
        return pdbqt_path
if __name__ == "__main__":
    args = parse()
    print(args)
    if args.smiles_path is not None:
        mol=[]
        smiles_list = []
        with open(args.smiles_path, "r") as f:
            for line in f:
                if line.strip() == "":
                    continue
                smiles_list.append(line.strip())
        #print(smiles_list)
        for i in range(len(smiles_list)):
            smiles = smiles_list[i]
            mol_name = f"molecule_{i}"
            try:
                mol.append(smiles_to_3D(smiles, mol_name))
            except Exception as e:
                print(f"Error converting smiles to 3D: {e}")
                continue
    elif args.smiles_str is not None:
        mol = [smiles_to_3D(args.smiles_str, "molecule_0")]
        
    os.makedirs("sdf", exist_ok=True)
    os.makedirs("pdbqt", exist_ok=True)
    for molecule in mol:
        Chem.MolToMolFile(molecule, f"sdf/{molecule.GetProp('_Name')}.sdf")
        pdbqt_path = sdf_to_pdbqt(f"sdf/{molecule.GetProp('_Name')}.sdf", f"pdbqt/{molecule.GetProp('_Name')}.pdbqt")
        print(pdbqt_path)
        