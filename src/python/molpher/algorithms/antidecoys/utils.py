import time

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate

from .settings import SIG_FAC

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds

def antifingerprint_from_paths(path, antifp_old=None):
    antifp_new = antifp_old
    for smiles in path:
        mol = Chem.MolFromSmiles(smiles)
        if not antifp_new:
            antifp_new = Generate.Gen2DFingerprint(mol, SIG_FAC)
        else:
            antifp_new = antifp_new | Generate.Gen2DFingerprint(mol, SIG_FAC)

    return antifp_new