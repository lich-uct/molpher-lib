import time

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds

class BidirectionalPathFinder:

    def __init__(self, source, target):
        options = {
            'fingerprint' : 'ATOM_PAIRS'
        }
        self.tree = ETree(source=source, target=target)
        self.ITERATION = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

    @property
    def path(self):
        current = self.tree.fetchMol(self.tree.params['target'])
        path = list()
        path.append(current.getSMILES())
        while current != '':
            current = current.getParentSMILES()
            if current:
                current = self.tree.fetchMol(current)
                path.append(current.getSMILES())

        path.reverse()
        return path

    def __call__(self):
        counter = 0
        connecting_molecule = None
        while not self.tree.path_found:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                print('Execution time ({0}):'.format(type(oper).__name__))

                time_elapsed = timeit(lambda : self.tree.runOperation(oper))
                print('\telapsed time: {0}'.format(time_elapsed))

def main():
    milliseconds_now = 1000 * time.clock()
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    pathfinder = BidirectionalPathFinder(cocaine, procaine)
    pathfinder()
    print(pathfinder.path)
    print('Total Execution Time: {0}'.format(1000 * time.clock() - milliseconds_now))

if __name__ == "__main__":
    exit(main())