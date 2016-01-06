from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

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

class BidirectionalPathFinder:

    def __init__(self, source, target):
        self.source_target = ETree(source=source, target=target)
        self.target_source = ETree(source=target, target=source)
        self.source_target.params = {
            'fingerprint' : 'ATOM_PAIRS'
        }
        self.target_source.params = {
            'fingerprint' : 'ATOM_PAIRS'
        }
        self.source_target_min = FindClosest()
        self.target_source_min = FindClosest()
        self.ITERATION = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

    def __call__(self):
        counter = 0
        while True:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                self.source_target.runOperation(oper)
                self.target_source.runOperation(oper)

            self.source_target.traverse(self.source_target_min)
            self.target_source.traverse(self.target_source_min)

            print('Current Targets:')
            print('source to target:', self.source_target.params['target'])
            print('target to source:', self.target_source.params['target'])

            print('Current Minima:')
            print('source to target:', self.source_target_min.closest.getSMILES(), self.source_target_min.closest.getDistToTarget())
            print('target to source:', self.target_source_min.closest.getSMILES(), self.target_source_min.closest.getDistToTarget())

            self.source_target.params = {
                'target' : self.target_source_min.closest.getSMILES()
            }
            self.target_source.params = {
                'target' : self.source_target_min.closest.getSMILES()
            }

            print('New Targets:')
            print('source to target:', self.source_target.params['target'])
            print('target to source:', self.target_source.params['target'])

            if self.source_target.path_found:
                print('Path Found in tree going from source to target:')
                connecting_molecule = self.source_target.params['target']
                print('Connecting molecule:', connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                break
            if self.target_source.path_found:
                print('Path Found in tree going from target to source:')
                connecting_molecule = self.target_source.params['target']
                print('Connecting molecule:', connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
                break

cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

pathfinder = BidirectionalPathFinder(cocaine, procaine)
pathfinder()
