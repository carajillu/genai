import argparse
import pandas as pd
import mdtraj
from dockflow.dock.prepare.ligands import smiles_to_pdbqt, smiles_to_mdtraj, sdf_to_pdbqt
from dockflow.dock.prepare.receptors import gen_binding_site_pdb, gen_receptor_pdbqt
from dockflow.dock.vina import gen_vina_config, dock
from dockflow.analyze.contact import pdbqt_to_mdtraj_trajectory, compute_contact_surface_area_single_pose
from dockflow.analyze.vina_scores import scores_from_pdbqt
from dockflow.analyze.rmsd import compute_rmsds_against_ref
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Run docking")
    parser.add_argument("--receptor", type=str, required=True, help="Path to the receptor PDB file")
    parser.add_argument("--library", type=str, default=None, help="Path to the ligand library of smiles files")
    parser.add_argument("--bsite_selection", type=str,required=True, help="region of the receptor we are gonna dock in")
    parser.add_argument("--ref_ligand", type=str, default=None, help="If relevant, reference ligand for rmsd calculation")
    args = parser.parse_args()
    if args.library is None and args.ref_ligand is None:
        raise ValueError("Either library or ref_ligand must be provided")
    return args

def main():
    args = parse_args()
    if args.library is not None:
        library = pd.read_csv(args.library)
    else:
        library = None
    root_dir = os.getcwd()
    receptor=mdtraj.load(args.receptor)
    if args.ref_ligand is not None:
        ref_ligand=mdtraj.load(args.ref_ligand.split(".")[0]+".pdb") #remove that and fix
    else:
        ref_ligand=None
    binding_site=gen_binding_site_pdb(receptor, args.bsite_selection)
    gen_receptor_pdbqt(args.receptor, "binding_site.pdb", "receptor")

    if args.ref_ligand is not None:
        os.makedirs("ref_ligand", exist_ok=True)
        os.chdir("ref_ligand")
        gen_vina_config(box_filename="../receptor.box.txt",out_path="vina_config.txt")
        cmd=f"mk_prepare_ligand.py -i {root_dir}/{args.ref_ligand} -o ref_ligand.pdbqt"
        os.system(cmd)
        dock(vina_exec="vina",receptor_pdbqt="../receptor.pdbqt",ligand_pdbqt=f"ref_ligand.pdbqt",config="vina_config.txt")
        csa=[]
        poses_mdtraj=pdbqt_to_mdtraj_trajectory(pdbqt_path="vina_out.pdbqt")
        for frame in poses_mdtraj:
            csa_frame=compute_contact_surface_area_single_pose(bsite_traj=binding_site,ligand_traj=frame)
            csa.append(csa_frame)
        scores=scores_from_pdbqt("vina_out.pdbqt")
        results=pd.DataFrame()
        results["ref_ligand"]=[args.ref_ligand]*len(scores)
        results["vina_rank"]=range(len(scores))
        results["vina_score"]=scores
        results["CSA"]=csa
        rmsd=compute_rmsds_against_ref(ref_pdbqt_path=f"ref_ligand.pdbqt",out_pdbqt_path="vina_out.pdbqt")
        results["rmsd"]=rmsd
        results.to_csv(f"ref_ligand.csv",index=False)
        os.chdir(root_dir)
    sys.exit()

    for index, row in library.iterrows():
        #prep
        smiles = row["smiles"]
        name = row["name"]
        os.makedirs(name, exist_ok=True)
        os.chdir(name)
        smiles_to_pdbqt(smiles=smiles, output_path=f"{name}.pdbqt",name=name)
        #dock
        gen_vina_config(box_filename="../receptor.box.txt",out_path="vina_config.txt")
        dock(vina_exec="vina",receptor_pdbqt="../receptor.pdbqt",ligand_pdbqt=f"{name}.pdbqt",config="vina_config.txt")
        #analysis
        csa=[]
        poses_mdtraj=pdbqt_to_mdtraj_trajectory(pdbqt_path="vina_out.pdbqt")
        for frame in poses_mdtraj:
            csa_frame=compute_contact_surface_area_single_pose(bsite_traj=binding_site,ligand_traj=frame)
            csa.append(csa_frame)
        scores=scores_from_pdbqt("vina_out.pdbqt")
        results=pd.DataFrame()
        results["name"]=[name]*len(scores)
        results["smiles"]=[smiles]*len(scores)
        results["vina_rank"]=range(len(scores))
        results["vina_score"]=scores
        results["CSA"]=csa
        if ref_ligand is not None:
            rmsd=compute_rmsds_against_ref(pdb_path=f"{root_dir}/{args.ref_ligand}",pdbqt_path="vina_out.pdbqt")
            results["rmsd"]=rmsd
        results.to_csv(f"{name}.csv",index=False)
        os.chdir(root_dir)

if __name__ == "__main__":
    main()