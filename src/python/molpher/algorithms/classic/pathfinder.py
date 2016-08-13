from molpher.algorithms.commons import find_path
from molpher.algorithms.operations import FindClosest
from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

class ClassicPathFinder:

    def __init__(self, settings):
        self.settings = settings
        self.tree = ETree.create(source=self.settings.source, target=self.settings.target)
        if self.settings.tree_params:
            self.tree.params = self.settings.tree_params
        self.tree.thread_count = self.settings.max_threads
        self.ITERATION = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper(settings.verbose)
            , ExtendTreeOper()
            , PruneTreeOper()
        ]
        self.find_closest = FindClosest()
        self.path = None

    def __call__(self):
        counter = 0
        while not self.tree.path_found:
            if counter > self.settings.max_iters:
                break
            counter+=1
            print('Iteration {0}'.format(counter))
            for oper in self.ITERATION:
                self.tree.runOperation(oper)
            self.tree.traverse(self.find_closest)
            print('Closest to target: {0} {1}'.format(self.find_closest.closest.smiles, self.find_closest.closest.dist_to_target))

        self.path = find_path(self.tree, self.tree.params['target'])
        print('Path found:', self.path)
        return self.path