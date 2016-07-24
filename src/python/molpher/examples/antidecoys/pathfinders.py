from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import ExtendTreeOper
from molpher.core.operations import FilterMorphsOper
from molpher.core.operations import GenerateMorphsOper
from molpher.core.operations import PruneTreeOper
from molpher.core.operations import SortMorphsOper

from .utils import timeit
from .settings import *
from .custom_opers import FindClosest

class BidirectionalPathFinder:

    def __init__(self, source, target, verbose=True, antidecoys_filter=None):
        self.antidecoys_filter = antidecoys_filter
        self.verbose = verbose

        self.source_target = ETree.create(source=source, target=target)
        self.source_target.thread_count = MAX_THREADS

        self.target_source = ETree.create(source=target, target=source)
        self.target_source.thread_count = MAX_THREADS

        self.source_target_min = FindClosest()
        self.target_source_min = FindClosest()

        print("Tree Parameters:")
        print('\tsource -> target: {0}'.format(self.source_target.params))
        print('\ttarget -> source: {0}'.format(self.target_source.params))

        self._iteration = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , self.antidecoys_filter
            , ExtendTreeOper()
            , PruneTreeOper()
        ]
        self._iteration = [x for x in self._iteration if x]
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
        search_failed = False
        antidecoys_off = False
        while True:
            counter+=1
            if counter > MAX_ITERS_PER_PATH:
                search_failed = True
                break
            if counter < WAIT_FOR_ANTIDECOYS + 1:
                antidecoys_off = True
            else:
                antidecoys_off = False
            print('Iteration {0}:'.format(counter))
            for oper in self._iteration:
                if oper == self.antidecoys_filter and antidecoys_off:
                    continue
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
                if issubclass(oper.__class__, GenerateMorphsOper):
                    print("Genereated morphs:")
                    print('\tsource -> target: {0}'.format(len(self.source_target.candidates)))
                    print('\ttarget -> source: {0}'.format(len(self.target_source.candidates)))

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

            source_target_min_dist = self.source_target_min.closest.getDistToTarget()
            target_source_min_dist = self.target_source_min.closest.getDistToTarget()
            print('Current Minima:')
            print('\tsource to target:', self.source_target_min.closest.getSMILES(), source_target_min_dist)
            print('\ttarget to source:', self.target_source_min.closest.getSMILES(), target_source_min_dist)
            if min(source_target_min_dist, target_source_min_dist) < ANTIDECOYS_DISTANCE_SWITCH:
                antidecoys_off = True
                print("Antidecoys turned off.")

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

        if not search_failed:
            source_target_path = self._find_path(self.source_target, connecting_molecule)
            target_source_path = self._find_path(self.target_source, connecting_molecule)
            assert source_target_path.pop(-1) == connecting_molecule
            target_source_path.reverse()
            source_target_path.extend(target_source_path)
            self.path = source_target_path
            print('Path:', self.path)
            return self.path
        else:
            print('Search reached maximum number of iterations. Aborting...')
            return None