import time

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate

from .settings import SIG_FAC, COMMON_BITS_PERC_THRS

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

def evaluate_path(path, anti_fp):
    path_fp = antifingerprint_from_paths(path)
    common_fp = path_fp & anti_fp

    on_bits = path_fp.GetNumOnBits()
    common_bits = common_fp.GetNumOnBits()
    common_bits_perc = common_bits / on_bits
    assert common_bits_perc <= 1.0
    if common_bits_perc > COMMON_BITS_PERC_THRS:
        return False, common_bits_perc
    else:
        return True, common_bits_perc

def compute_antifp_scores(mols, anti_fp):
    scores = dict()
    for mol in mols:
        scores[mol] = evaluate_path([mol], anti_fp)[1]
    return scores