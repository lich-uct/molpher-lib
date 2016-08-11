import time

import multiprocessing
from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds

def compute_anti_fp(mols_smiles, sig_fac, antifp_old=None):
    antifp_new = antifp_old
    for smiles in mols_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if not antifp_new:
            antifp_new = Generate.Gen2DFingerprint(mol, sig_fac)
        else:
            antifp_new = antifp_new | Generate.Gen2DFingerprint(mol, sig_fac)

    return antifp_new

def eval_morphs(mols, anti_fp, sig_fac):
    path_fp = compute_anti_fp(mols, sig_fac)
    common_fp = path_fp & anti_fp

    on_bits = path_fp.GetNumOnBits()
    common_bits = common_fp.GetNumOnBits()
    common_bits_perc = common_bits / on_bits if on_bits else 0.0
    assert common_bits_perc <= 1.0
    return common_bits_perc

def compute_antifp_score(data):
    smiles, anti_fp = data
    return smiles, eval_morphs([smiles], anti_fp, compute_antifp_score.sig_fac)

def compute_antifp_scores(mols, anti_fp, settings):
    data = [(x.smiles, anti_fp) for x in mols]

    compute_antifp_score.sig_fac = settings.signature_factory
    pool = multiprocessing.Pool(settings.max_threads)
    results = pool.map(func=compute_antifp_score, iterable=data)
    pool.close()
    pool.join()
    return dict(results)

def update_target(tree, target):
    if target != tree.params['source']:
        tree.params = {
            'target' : target
        }

def find_path(tree, end_mol):
    path = []
    current = tree.fetchMol(end_mol)
    path.append(current.getSMILES())
    while current != '':
        current = current.getParentSMILES()
        if current:
            current = tree.fetchMol(current)
            path.append(current.getSMILES())
    path.reverse()
    return path