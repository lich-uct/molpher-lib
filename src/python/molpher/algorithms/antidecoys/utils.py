import time
import pickle

from .settings import *
from .custom_opers import AntidecoysFilterMulti

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate

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

def run(source, target, antidecoys_filter=AntidecoysFilterMulti()):
    antifp_path = os.path.join(STORAGE_DIR, 'antifingerprint.pickle')
    paths_path = os.path.join(STORAGE_DIR, 'paths.pickle')
    from .pathfinders import BidirectionalPathFinder

    paths = []
    if os.path.exists(paths_path):
        pickled_paths = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_paths))
        pickled_paths.close()
    if os.path.exists(antifp_path):
        pickled_antifp = open(antifp_path, mode='rb')
        antidecoys_filter.antifingerprint = pickle.load(pickled_antifp)
        pickled_antifp.close()
    for i in range(PATHS_TO_FIND):
        # find a path
        pathfinder = BidirectionalPathFinder(source, target, verbose=False, antidecoys_filter=antidecoys_filter)
        exec_time = timeit(pathfinder)
        paths.append(pathfinder.path)
        print('Total Execution Time (search #{1}): {0}'.format(exec_time, i + 1))

        # compute and save new antifingerprint
        antidecoys_filter.antifingerprint = antifingerprint_from_paths(pathfinder.path, antidecoys_filter.antifingerprint)

        # pickle the results for future use
        pickled_antifp = open(antifp_path, mode='wb')
        pickle.dump(antidecoys_filter.antifingerprint, pickled_antifp)
        pickled_antifp.close()
        pickled_paths = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_paths)
        pickled_paths.close()