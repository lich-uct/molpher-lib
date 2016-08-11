from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

class ClassicPathFinder:

    def __init__(self, source, target):
        self.tree = ETree.create(source=source, target=target)
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
        path.append(current.smiles)
        while current != '':
            current = current.parent_smiles
            if current:
                current = self.tree.fetchMol(current)
                path.append(current.smiles)

        path.reverse()
        return path

    def __call__(self):
        counter = 0
        while not self.tree.path_found:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                self.tree.runOperation(oper)