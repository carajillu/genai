import MDAnalysis as mda
from MDAnalysis.analysis import rms
import os
from meeko import PDBQTMolecule

def compute_rmsds_against_ref(ref_pdbqt_path, out_pdbqt_path,n_poses=1000):
    # Load reference ligand (PDB) with MDTraj and remove Hs
    ref = mda.Universe(ref_pdbqt_path)
    # Load PDBQT file and extract all poses (MODEL blocks)
    poses = PDBQTMolecule.from_file(out_pdbqt_path)
    rmsd_list = []
    for pose in poses:
        pose.write_pdbqt_file("pose.pdbqt")
        pose = mda.Universe("pose.pdbqt")
        rmsd = rms.rmsd(pose.atoms.positions,ref.atoms.positions,center=False,superposition=False)/10 #convert to nm
        rmsd_list.append(rmsd)
        os.remove("pose.pdbqt")
    return rmsd_list

