import subprocess
import mdtraj
import sys

def gen_binding_site_pdb(receptor_path: str="receptor.pdb", selection: str="all"):
    """
    Select the binding site of a receptor
    """
    receptor= mdtraj.load(receptor_path)
    binding_site = receptor.atom_slice(receptor.topology.select(selection))
    binding_site_filename = f"binding_site.pdb"
    binding_site.save_pdb(binding_site_filename)
    return binding_site_filename

def gen_receptor_pdbqt(receptor_path: str, binding_site_pdb: str, output_basename: str = "receptor"):
    """
    Generate a pdbqt file for a receptor and a docking box
    TODO: use proper python functions from meeko
    """
    cmd_str=f"mk_prepare_receptor.py --read_pdb {receptor_path} -o {output_basename} -p -v --box_enveloping {binding_site_pdb} --padding 0 -a"
    print(f"Running command: {cmd_str}")
    cmd=cmd_str.split()
    subprocess.run(cmd)
    return 

