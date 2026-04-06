from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np

def internal_diversity(smiles, radius=2, nBits=2048):
    mols = [Chem.MolFromSmiles(s) for s in smiles if Chem.MolFromSmiles(s)]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits) for m in mols]

    sims = []
    for i in range(len(fps)):
        sims.extend(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:]))

    return 1-np.mean(sims)
