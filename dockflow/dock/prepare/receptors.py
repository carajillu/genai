import subprocess
import mdtraj

def gen_binding_site_pdb(receptor: mdtraj.Trajectory, selection: str):
    """
    Select the binding site of a receptor
    """
    binding_site = receptor.atom_slice(receptor.topology.select(selection))
    binding_site_filename = f"binding_site.pdb"
    binding_site.save_pdb(binding_site_filename)
    return binding_site

def gen_receptor_pdbqt(receptor_pdb: str, binding_site_pdb: str, output_path: str = "receptor"):
    """
    Generate a pdbqt file for a receptor and a docking box
    TODO: use proper python functions from meeko
    """
    cmd_str=f"mk_prepare_receptor.py --read_pdb {receptor_pdb} -o {output_path} -p -v --box_enveloping {binding_site_pdb} --padding 0 -a"
    cmd=cmd_str.split()
    subprocess.run(cmd)
    return

