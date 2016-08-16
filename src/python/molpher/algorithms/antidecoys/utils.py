import multiprocessing

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate


def compute_anti_fp(mols_smiles, sig_fac, antifp_old=None):
    """
    Computes an anti-fingerprint from the given molecules.

    It is possible to specify an existing anti-fingerprint
    for an update. In this case, the returned fingerprint
    will be the result of doing a bitwise :samp:`or` between
    the old fingerprint and the one generated from the
    supplied structures.

    :param mols_smiles: SMILES of the molecules to generate the anti-fingerprint from
    :param sig_fac: RDKit's signature factory used in the 2D pharmacophore fingerprint computation
    :param antifp_old: an old anti-fingerprint to update
    :return: new or updated anti-fingerprint
    """

    antifp_new = antifp_old
    for smiles in mols_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if not antifp_new:
            antifp_new = Generate.Gen2DFingerprint(mol, sig_fac)
        else:
            antifp_new = antifp_new | Generate.Gen2DFingerprint(mol, sig_fac)

    return antifp_new

def _eval_morphs(mols, anti_fp, sig_fac):
    path_fp = compute_anti_fp(mols, sig_fac)
    common_fp = path_fp & anti_fp

    on_bits = path_fp.GetNumOnBits()
    common_bits = common_fp.GetNumOnBits()
    common_bits_perc = common_bits / on_bits if on_bits else 0.0
    assert common_bits_perc <= 1.0
    return common_bits_perc

def _compute_antifp_score(data):
    smiles, anti_fp = data
    return smiles, _eval_morphs([smiles], anti_fp, _compute_antifp_score.sig_fac)

def compute_antifp_scores(mols, anti_fp, settings):
    """
    Computes (in parallel) the anti-fingerprint scores
    for the given molecules.

    :param mols: molecules as a list of SMILES strings
    :param anti_fp: the anti-fingerprint
    :param settings: `AntidecoysSettings`
    :return: a `dict` mapping SMILES of the given molecules to the anti-fingerprint scores
    """

    data = [(x.smiles, anti_fp) for x in mols]

    _compute_antifp_score.sig_fac = settings.signature_factory # RDKit's signature factory cannot be pickled; hence, this dirty hack
    pool = multiprocessing.Pool(settings.max_threads)
    results = pool.map(func=_compute_antifp_score, iterable=data)
    pool.close()
    pool.join()
    return dict(results)

def update_target(tree, target):
    """
    Change the tree's target molecule.

    :param tree: an `exploration tree`
    :type tree: :class:`~molpher.core.ExplorationTree.ExplorationTree`
    :param target: SMILES string of the new target
    :type target: `str`
    """

    if target != tree.params['source']:
        tree.params = {
            'target' : target
        }
