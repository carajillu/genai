import mdtraj
from vina import Vina
import meeko
import subprocess
import pandas as pd
import argparse
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", type=str, required=True)
    parser.add_argument("--ligand", type=str, required=True)
    parser.add_argument("--binding_site_selection", type=str, required=True)
    parser.add_argument("--vina_exec", type=str, default="vina")
    return parser.parse_args()

def gen_binding_site(receptor: mdtraj.Trajectory, selection: str):
    """
    Select the binding site of a receptor
    """
    binding_site = receptor.atom_slice(receptor.topology.select(selection))
    binding_site_filename = f"binding_site.pdb"
    binding_site.save_pdb(binding_site_filename)
    return

def gen_receptor_pdbqt(receptor_pdb: str, binding_site_pdb: str):
    """
    Generate a pdbqt file for a receptor and a docking box
    TODO: use proper python functions from meeko
    """
    cmd_str=f"mk_prepare_receptor.py --read_pdb {receptor_pdb} -o receptor -p -v --box_enveloping {binding_site_pdb} --padding 0 -a"
    cmd=cmd_str.split()
    subprocess.run(cmd)
    return

def gen_ligand_pdbqt(ligand_sdf: str):
    """
    Generate a pdbqt file for a ligand
    """
    cmd_str=f"mk_prepare_ligand.py -i {ligand_sdf} -o ligand.pdbqt"
    cmd=cmd_str.split()
    subprocess.run(cmd)

def run_docking(receptor_pdbqt: str, ligand_pdbqt: str, box_filename: str, vina_exec: str = "vina"):
    """
    Run docking
    """
    box_df=pd.read_csv(box_filename, sep=" ", header=None, names=["field","equal_sign","value"]) #column names are "field","equal_sign", and "value"
    box_center=list(box_df.value[0:3])
    box_size=list(box_df.value[3:6])
    
    with open("vina_config.txt", "w") as f:
        f.write(f"receptor = {receptor_pdbqt}\n")
        f.write(f"ligand = {ligand_pdbqt}\n")
        f.write(f"center_x = {box_center[0]}\n")
        f.write(f"center_y = {box_center[1]}\n")
        f.write(f"center_z = {box_center[2]}\n")
        f.write(f"size_x = {box_size[0]}\n")
        f.write(f"size_y = {box_size[1]}\n")
        f.write(f"size_z = {box_size[2]}\n")
        f.write(f"spacing = 0.375\n") 
        f.write(f"exhaustiveness = 32\n")
        f.write(f"num_modes = 100\n")
        f.write(f"out = vina_out.pdbqt\n")
        f.write(f"energy_range = 99999.0\n")
    cmd_str=f"{vina_exec} --config vina_config.txt"
    cmd=cmd_str.split()
    subprocess.run(cmd)
       
    return
    '''
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    box_df=pd.read_csv(box_filename, sep=" ", header=None)[2]
    v.compute_vina_maps(center=box_df[0:3], box_size=box_df[3:6])
    v.dock(exhaustiveness=32, n_poses=100)
    v.write_poses('vina_out.pdbqt', n_poses=100, overwrite=True)
    return
    '''


if __name__ == "__main__":
    args = parse_args()
    receptor = mdtraj.load_pdb(args.receptor)
    gen_binding_site(receptor, args.binding_site_selection)
    gen_receptor_pdbqt(args.receptor, "binding_site.pdb")
    gen_ligand_pdbqt(args.ligand)
    run_docking("receptor.pdbqt", "ligand.pdbqt", "receptor.box.txt", args.vina_exec)
