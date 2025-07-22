import os
import sys
import subprocess
import argparse
import pandas as pd
import numpy as np
import glob
import mdtraj
from dockflow.dock.prepare.ligands import prepare_library_from_smiles
from dockflow.dock.prepare.receptors import gen_binding_site_pdb, gen_receptor_pdbqt
from dockflow.dock.vina import gen_vina_config, batch_dock
from dockflow.analyze.contact import pdbqt_to_mdtraj_trajectory, compute_contact_surface_area_single_pose
from dockflow.analyze.vina_scores import scores_from_pdbqt
from dockflow.analyze.rmsd import compute_rmsds_against_ref
from meeko import PDBQTMolecule
def parse_args():
    parser = argparse.ArgumentParser(description="Run parallel docking")
    parser.add_argument("--receptor", type=str, help="Path to the receptor PDB file if only one")
    parser.add_argument("--receptor_dir", type=str, default=None, help="Path to the directory containing receptor PDB files")
    parser.add_argument("--library_smi", type=str, default=None, help="Path to the ligand library of smiles files")
    parser.add_argument("--bsite_selection", type=str, required=True, help="Region of the receptor we are going to dock in")
    parser.add_argument("--ref_ligand", type=str, default=None, help="If relevant, reference ligand for RMSD calculation")
    parser.add_argument("--score_threshold", type=float, default=-10.0, help="Score threshold for filtering docking results")
    parser.add_argument("--results", type=str, default="results.csv", help="Path to save the ligand candidates that passed")
    args = parser.parse_args()
    if args.library_smi is None and args.ref_ligand is None:
        raise ValueError("Either --library_smi or --ref_ligand must be provided")
    if args.receptor is None and args.receptor_dir is None:
        raise ValueError("Either --receptor or --receptor_dir must be provided")

    return args

if __name__=="__main__":
    
    args=parse_args()

    # Read SMILES file and prepare library
    smiles_df=pd.read_csv(args.library_smi)
    assert 'smiles' in smiles_df.columns or 'SMILES' in smiles_df.columns, \
        "SMILES file must contain a column named 'smiles' or 'SMILES'."
    smiles_list = smiles_df['smiles'] if 'smiles' in smiles_df.columns else smiles_df['SMILES']

    pdbqt_library_path=prepare_library_from_smiles(smiles_list, output_dir="pdbqt_library")
    best_scores=[np.inf]*len(smiles_list) # Initialize best scores for each ligand
    receptor_conformations=[None]*len(smiles_list) # receptor pdb that docks the best ligand pose

    # get receptor paths
    if args.receptor_dir is not None:
        receptor_paths = glob.glob(os.path.join(args.receptor_dir, "*.pdb"))
    elif args.receptor is not None:
        receptor_paths = [args.receptor]
    
    root_dir = os.getcwd()
    for receptor in receptor_paths:
        # Prepare receptor path and name
        receptor_path = os.path.abspath(receptor)
        receptor_name = receptor.split("/")[-1].split(".")[0]
        os.makedirs(receptor_name, exist_ok=True)
        os.chdir(receptor_name)

        # Prepare binding site PDB and box
        binding_site = gen_binding_site_pdb(receptor_path=receptor_path, selection=args.bsite_selection)
        
        # Prepare receptor PDBQT
        gen_receptor_pdbqt(receptor_path=receptor_path,binding_site_pdb=binding_site, output_basename=receptor_name)
        
        
        
        # dock the ligands in each receptor
        batch_dock(vina_exec="~/github/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1", receptor_pdbqt=f"{receptor_name}.pdbqt", 
                              ligand_directory=pdbqt_library_path, config=f"{receptor_name}.box.txt",
                              thread=8000,opencl_binary_path="/home/joan/github/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1",output_directory="./results")
        
        #pose_bust(vina_output_path) # run posebusters on all poses, discard the ones that don't pass

        for i in range(len(smiles_list)):
            ligand_out= f"results/ligand_{i}_out.pdbqt"
            poses=PDBQTMolecule.from_file(ligand_out)
            score=poses.score # that returns the score from the first pose, which is the best one
            if score < best_scores[i]:
                best_scores[i] = score
                receptor_conformations[i] = receptor_name
        os.chdir(root_dir)
    
    
    smiles_df['best_score'] = best_scores
    smiles_df['receptor'] = receptor_conformations
    smiles_df.to_csv(args.results, index=False)