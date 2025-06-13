import mdtraj as md
from mdtraj.core.topology import Topology

def renumber_and_split_chains(pdb_path, output_path):
    traj = md.load_pdb(pdb_path)
    old_top = traj.topology

    # Define offset based on given mappings
    chain_a_offset = -23  # all residues in 52–444
    chain_b_offset = 6    # all residues in 1–31

    # Separate residues
    residues = list(old_top.residues)
    chain_a_residues = [r for r in residues if 52 <= r.resSeq <= 444]
    chain_b_residues = [r for r in residues if 1 <= r.resSeq <= 31 and r.resSeq < 52]

    # Create new topology
    new_top = Topology()
    atom_mapping = {}
    new_chain_a = new_top.add_chain('A')
    new_chain_b = new_top.add_chain('B')
    atom_indices_order = []

    def copy_residues_with_offset(res_list, chain, offset):
        indices = []
        for res in res_list:
            new_resSeq = res.resSeq + offset
            new_res = new_top.add_residue(res.name, chain, resSeq=new_resSeq)
            for atom in res.atoms:
                new_atom = new_top.add_atom(atom.name, atom.element, new_res)
                atom_mapping[atom.index] = new_atom
                indices.append(atom.index)
        return indices

    # Copy and renumber
    atom_indices_a = copy_residues_with_offset(chain_a_residues, new_chain_a, chain_a_offset)
    atom_indices_b = copy_residues_with_offset(chain_b_residues, new_chain_b, chain_b_offset)

    # Combine indices and reorder coordinates
    new_atom_order = atom_indices_a + atom_indices_b
    new_coords = traj.xyz[:, new_atom_order, :]
    new_traj = md.Trajectory(xyz=new_coords, topology=new_top)

    # Save new PDB
    new_traj.save_pdb(output_path)

def find_nonstandard_residues(traj):
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
        'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO',
        'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

    nonstandard = set()

    for res in traj.topology.residues:
        if res.name.upper() not in standard_aa:
            nonstandard.add(res.name)

    return sorted(nonstandard)
