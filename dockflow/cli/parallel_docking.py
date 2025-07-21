import os
import sys
import subprocess
import argparse
import pandas as pd
import mdtraj
from dockflow.dock.prepare.ligands import prepare_library_from_smiles
from dockflow.dock.prepare.receptors import gen_binding_site_pdb, gen_receptor_pdbqt
from dockflow.dock.vina import gen_vina_config, batch_dock
from dockflow.analyze.contact import pdbqt_to_mdtraj_trajectory, compute_contact_surface_area_single_pose
from dockflow.analyze.vina_scores import scores_from_pdbqt
from dockflow.analyze.rmsd import compute_rmsds_against_ref
def parse_args():
    parser = argparse.ArgumentParser(description="Run parallel docking")
    parser.add_argument("--receptor",nargs="+", type=str, required=True, help="Path to the receptor PDB file")
    parser.add_argument("--library_smi", type=str, default=None, help="Path to the ligand library of smiles files")
    parser.add_argument("--bsite_selection", type=str, required=True, help="Region of the receptor we are going to dock in")
    parser.add_argument("--ref_ligand", type=str, default=None, help="If relevant, reference ligand for RMSD calculation")
    parser.add_argument("--score_threshold", type=float, default=-10.0, help="Score threshold for filtering docking results")
    parser.add_argument("--results", type=str, default="results.csv", help="Path to save the ligand candidates that passed")
    args = parser.parse_args()
    if args.library_smi is None and args.ref_ligand is None:
        raise ValueError("Either --library_smi or --ref_ligand must be provided")
    return args

if __name__=="__main__":
    
    args=parse_args()
    pdbqt_library_path=prepare_library_from_smiles(args.library_smi, output_dir="pdbqt_library")
    root_dir = os.getcwd()
    for receptor in args.receptor:
        receptor_path = os.path.abspath(receptor)
        receptor_name = receptor.split("/")[-1].split(".")[0]
        os.makedirs(receptor_name, exist_ok=True)
        os.chdir(receptor_name)

        # Prepare binding site PDB and box
        binding_site = gen_binding_site_pdb(receptor_path=receptor_path, selection=args.bsite_selection)
        
        # Prepare receptor PDBQT
        gen_receptor_pdbqt(receptor_path=receptor_path,binding_site_pdb=binding_site, output_basename=receptor_name)
        
        
        
        # dock the ligands in each receptor
        vina_output_path=batch_dock(vina_exec="~/github/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1", receptor_pdbqt=f"{receptor_name}.pdbqt", 
                              ligand_directory=pdbqt_library_path, config=f"{receptor_name}.box.txt",
                              thread=8000,opencl_binary_path="/home/joan/github/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1",output_directory="./results")
        
        pose_bust(vina_output_path) # run posebusters on all poses, discard the ones that don't pass

        results=get_results(vina_output_path, score_threshold, results) # get the results of each docking and append them to a DataFrame

        os.chdir(root_dir)
    sys.exit()
    best_results=consolidate_results() # remove redundant smiles, keep the best scoring ones, discard compounds above similarity threshold (?)
    best_results.to_csv("best_results.csv", index=False)  # Save the consolidated results
