import os
import pickle
import time

from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures

from rdkit import Chem
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import *
from molpher.core.selectors import *
from molpher import random

#random.set_random_seed(42)

# dir for stored data
STORAGE_DIR = os.path.abspath('data')
if not os.path.exists(STORAGE_DIR):
    os.mkdir(STORAGE_DIR)

# number of threads to use
THREADS = 2

# important stuff for the pharmacophore fingerprints
FDEF_FILE = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')  # get basic feature definitions
FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(FDEF_FILE)  # make feature factory
SIG_FAC = SigFactory(FEATURE_FACTORY, minPointCount=2, maxPointCount=3, trianglePruneBins=False)  # make signature factory
SIG_FAC.SetBins([(0, 2), (2, 5), (5, 8)])
SIG_FAC.Init()

# threshold for the bits in common percentage
COMMON_BITS_PERC_THRS = 0.5

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds

class FindClosest:

    def __init__(self):
        self.closest = None

    def __call__(self, morph):
        if not self.closest:
            self.closest = morph.copy()
            return
        current_dist = self.closest.getDistToTarget()
        morph_dist = morph.getDistToTarget()
        if morph_dist < current_dist:
            self.closest = morph.copy()

class AntidecoysFilter(TreeOperation):

    def __init__(self):
        super(AntidecoysFilter, self).__init__()
        self.antifingerprint = None

    def __call__(self):
        if self.antifingerprint:
            tree = self.tree
            mask = list(tree.candidates_mask)
            candidates = [Chem.MolFromSmiles(mol.smiles) if mask[idx] else None for idx, mol in enumerate(tree.candidates)]
            candidates_fps = [Generate.Gen2DFingerprint(mol, SIG_FAC) if mol else None for mol in candidates]

            counter = 0
            for idx, fp in enumerate(candidates_fps):
                if fp:
                    candidate_bits = fp.GetNumOnBits()
                    common_fp = self.antifingerprint & fp
                    common_bits = common_fp.GetNumOnBits()
                    common_bits_perc = common_bits / candidate_bits
                    if common_bits_perc > COMMON_BITS_PERC_THRS:
                        mask[idx] = False
                        counter += 1
            print('Anti-decoys filter eliminated {0} molecules...'.format(counter))
            tree.candidates_mask = mask


class BidirectionalPathFinder:

    def __init__(self, source, target, verbose=True, antidecoys_filter=AntidecoysFilter()):
        self.antidecoys_filter = antidecoys_filter
        self.verbose = verbose

        self.source_target = ETree.create(source=source, target=target)
        self.source_target.thread_count = THREADS

        self.target_source = ETree.create(source=target, target=source)
        self.target_source.thread_count = THREADS

        self.source_target_min = FindClosest()
        self.target_source_min = FindClosest()

        self.ITERATION = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , self.antidecoys_filter
            , ExtendTreeOper()
            , PruneTreeOper()
        ]
        self.path = []

    def _find_path(self, tree, connecting_mol):
        path = []
        current = tree.fetchMol(connecting_mol)
        path.append(current.getSMILES())
        while current != '':
            current = current.getParentSMILES()
            if current:
                current = tree.fetchMol(current)
                path.append(current.getSMILES())
        path.reverse()
        return path

    def __call__(self):
        counter = 0
        connecting_molecule = None
        while True:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                if self.verbose:
                    print('Execution times ({0}):'.format(type(oper).__name__))

                    source_target_time = timeit(lambda : self.source_target.runOperation(oper))
                    print('\tsource -> target: {0}'.format(source_target_time))
                    target_source_time = timeit(lambda : self.target_source.runOperation(oper))
                    print('\ttarget -> source: {0}'.format(target_source_time))

                    print('\ttotal time: {0}'.format(source_target_time + target_source_time))
                else:
                    self.source_target.runOperation(oper)
                    self.target_source.runOperation(oper)

            if self.verbose:
                print('Traversal times:')

                source_target_time = timeit(lambda : self.source_target.traverse(self.source_target_min))
                print('\tsource -> target: {0}'.format(source_target_time))
                target_source_time = timeit(lambda : self.target_source.traverse(self.target_source_min))
                print('\ttarget -> source: {0}'.format(target_source_time))

                print('\ttotal execution time: {0}'.format(source_target_time + target_source_time))

                print('Current Targets:')
                print('\tsource to target:', self.source_target.params['target'])
                print('\ttarget to source:', self.target_source.params['target'])
            else:
                self.source_target.traverse(self.source_target_min)
                self.target_source.traverse(self.target_source_min)

            print('Current Minima:')
            print('\tsource to target:', self.source_target_min.closest.getSMILES(), self.source_target_min.closest.getDistToTarget())
            print('\ttarget to source:', self.target_source_min.closest.getSMILES(), self.target_source_min.closest.getDistToTarget())

            self.source_target.params = {
                'target' : self.target_source_min.closest.getSMILES()
            }
            self.target_source.params = {
                'target' : self.source_target_min.closest.getSMILES()
            }

            if self.verbose:
                print('New Targets:')
                print('\tsource to target:', self.source_target.params['target'])
                print('\ttarget to source:', self.target_source.params['target'])

            if self.source_target.path_found:
                connecting_molecule = self.source_target.params['target']
                if self.verbose:
                    print('Path Found in tree going from source to target')
                    print('Connecting molecule:', connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                break
            if self.target_source.path_found:
                connecting_molecule = self.target_source.params['target']
                if self.verbose:
                    print('Path Found in tree going from target to source')
                    print('Connecting molecule:', connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
                break

        source_target_path = self._find_path(self.source_target, connecting_molecule)
        target_source_path = self._find_path(self.target_source, connecting_molecule)
        assert source_target_path.pop(-1) == connecting_molecule
        target_source_path.reverse()
        source_target_path.extend(target_source_path)
        self.path = source_target_path
        print('Path:', self.path)
        return self.path


def antifingerprint_from_paths(path, antifp_old=None):
    antifp_new = antifp_old
    for smiles in path:
        mol = Chem.MolFromSmiles(smiles)
        if not antifp_new:
            antifp_new = Generate.Gen2DFingerprint(mol, SIG_FAC)
        else:
            antifp_new = antifp_new | Generate.Gen2DFingerprint(mol, SIG_FAC)

    return antifp_new

def main():
    milliseconds_now = 1000 * time.clock()
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    antifp_path = os.path.join(STORAGE_DIR, 'antifingerprint.pickle')
    paths_path = os.path.join(STORAGE_DIR, 'paths.pickle')

    paths = []
    if os.path.exists(paths_path):
        pickled_paths = open(paths_path, mode='rb')
        paths.extend(pickle.load(pickled_paths))
        pickled_paths.close()
    antidecoys_filter = AntidecoysFilter()
    if os.path.exists(antifp_path):
        pickled_antifp = open(antifp_path, mode='rb')
        antidecoys_filter.antifingerprint = pickle.load(pickled_antifp)
        pickled_antifp.close()
    for i in range(5):
        # find a path
        pathfinder = BidirectionalPathFinder(cocaine, procaine, verbose=False, antidecoys_filter=antidecoys_filter)
        pathfinder()
        paths.append(pathfinder.path)
        print('Total Execution Time (search #{1}): {0}'.format(1000 * time.clock() - milliseconds_now, i + 1))

        # compute and save new antifingerprint
        antidecoys_filter.antifingerprint = antifingerprint_from_paths(pathfinder.path, antidecoys_filter.antifingerprint)

        # pickle the results for future use
        pickled_antifp = open(antifp_path, mode='wb')
        pickle.dump(antidecoys_filter.antifingerprint, pickled_antifp)
        pickled_antifp.close()
        pickled_paths = open(paths_path, mode='wb')
        pickle.dump(paths, pickled_paths)
        pickled_paths.close()


if __name__ == "__main__":
    exit(main())