import subprocess
import pandas as pd
import os

def gen_vina_config(box_filename: str = "receptor.box.txt", out_path: str = "vina_config.txt"):
    box_df=pd.read_csv(box_filename, sep=" ", header=None, names=["field","equal_sign","value"]) #column names are "field","equal_sign", and "value"
    box_center=list(box_df.value[0:3])
    box_size=list(box_df.value[3:6])
    with open("vina_config.txt", "w") as f:
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
    return

def dock(vina_exec: str="vina", receptor_pdbqt:str = "receptor.pdbqt", ligand_pdbqt: str = "ligand.pdbqt", config="vina_config.txt", **kwargs):
    """
    Run Vina docking with the specified parameters.
    """
    cmd=f"{vina_exec} --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --config {config}"
    if kwargs:
        for key, value in kwargs.items():
            cmd += f" --{key} {value}"
    print(cmd)
    os.system(cmd)
    #subprocess.run(cmd.split(), shell=True)

def batch_dock(vina_exec: str="vina", receptor_pdbqt:str = "receptor.pdbqt", ligand_directory: str = "ligands/", config="vina_config.txt", **kwargs):
    cmd=f"{vina_exec} --receptor {receptor_pdbqt} --ligand_directory {ligand_directory} --config {config}"
    if kwargs:
        for key, value in kwargs.items():
            cmd += f" --{key} {value}"
    print(cmd)
    os.system(cmd)