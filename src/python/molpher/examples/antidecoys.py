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

random.set_random_seed(42)

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
            # TODO: implement
            pass

class BidirectionalPathFinder:

    def __init__(self, source, target, verbose=True, antidecoys_filter=AntidecoysFilter()):
        self.antidecoys_filter = antidecoys_filter
        self.verbose = verbose
        options = {
            'fingerprint' : FP_ATOM_PAIRS
        }
        self.source_target = ETree.create(source=source, target=target)
        self.target_source = ETree.create(source=target, target=source)
        self.source_target.params = options
        self.target_source.params = options
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
                print('Execution times ({0}):'.format(type(oper).__name__)) if self.verbose else None

                source_target_time = timeit(lambda : self.source_target.runOperation(oper))
                print('\tsource -> target: {0}'.format(source_target_time)) if self.verbose else None
                target_source_time = timeit(lambda : self.target_source.runOperation(oper))
                print('\ttarget -> source: {0}'.format(target_source_time)) if self.verbose else None

                print('\ttotal time: {0}'.format(source_target_time + target_source_time)) if self.verbose else None

            print('Traversal times:') if self.verbose else None

            source_target_time = timeit(lambda : self.source_target.traverse(self.source_target_min))
            print('\tsource -> target: {0}'.format(source_target_time)) if self.verbose else None
            target_source_time = timeit(lambda : self.target_source.traverse(self.target_source_min))
            print('\ttarget -> source: {0}'.format(target_source_time)) if self.verbose else None

            print('\ttotal execution time: {0}'.format(source_target_time + target_source_time))

            print('Current Targets:') if self.verbose else None
            print('\tsource to target:', self.source_target.params['target']) if self.verbose else None
            print('\ttarget to source:', self.target_source.params['target']) if self.verbose else None

            print('Current Minima:')
            print('\tsource to target:', self.source_target_min.closest.getSMILES(), self.source_target_min.closest.getDistToTarget())
            print('\ttarget to source:', self.target_source_min.closest.getSMILES(), self.target_source_min.closest.getDistToTarget())

            self.source_target.params = {
                'target' : self.target_source_min.closest.getSMILES()
            }
            self.target_source.params = {
                'target' : self.source_target_min.closest.getSMILES()
            }

            print('New Targets:') if self.verbose else None
            print('\tsource to target:', self.source_target.params['target']) if self.verbose else None
            print('\ttarget to source:', self.target_source.params['target']) if self.verbose else None

            if self.source_target.path_found:
                print('Path Found in tree going from source to target') if self.verbose else None
                connecting_molecule = self.source_target.params['target']
                print('Connecting molecule:', connecting_molecule) if self.verbose else None
                assert self.source_target.hasMol(connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                break
            if self.target_source.path_found:
                print('Path Found in tree going from target to source') if self.verbose else None
                connecting_molecule = self.target_source.params['target']
                print('Connecting molecule:', connecting_molecule) if self.verbose else None
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


def antifingerprint_from_paths(path, antifingerprint=None):
    fdef_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')  # get basic feature definitions
    feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_name)  # make feature factory
    sign_fac = SigFactory(feature_factory, minPointCount=2, maxPointCount=3)
    sign_fac.SetBins([(0,2),(2,5),(5,8)])
    sign_fac.Init()

    antifingerprint = antifingerprint
    for smiles in path:
        mol = Chem.MolFromSmiles(smiles)
        if not antifingerprint:
            antifingerprint = Generate.Gen2DFingerprint(mol, sign_fac)
        else:
            antifingerprint = antifingerprint | Generate.Gen2DFingerprint(mol, sign_fac)

    return antifingerprint

def main():
    milliseconds_now = 1000 * time.clock()
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    paths = []
    antidecoys_filter = AntidecoysFilter()
    for i in range(10):
        # find a path
        pathfinder = BidirectionalPathFinder(cocaine, procaine, verbose=False, antidecoys_filter=antidecoys_filter)
        paths.append(pathfinder.path)
        pickle.dump(paths, open('paths.pickle'.format(i)), mode='bw')
        print('Total Execution Time (search #{1}): {0}'.format(1000 * time.clock() - milliseconds_now, i + 1))

        # compute and save new antifingerprint
        antidecoys_filter.antifingerprint = antifingerprint_from_paths(pathfinder.path, antidecoys_filter.antifingerprint)


if __name__ == "__main__":
    exit(main())