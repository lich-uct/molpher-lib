from molpher.algorithms.functions import find_path, timeit
from molpher.algorithms.operations import FindClosest
from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import *

class BidirectionalPathFinder:
    """
    :param settings: search settings
    :type settings: `Settings`

    Implements a search where two trees are used to explore against each other.
    After every iteration the targets of the trees are updated so that the closest molecule
    from one tree becomes a target in the other.

    This search might be quicker and the morphs on the resulting paths are more balanced
    in terms of similarity to both source and the target.
    """

    def __init__(self, settings, iter_callback=None):
        self.verbose = settings.verbose
        """`True` if verbose output was requested"""

        self.iter_callback = iter_callback
        """callback to call after each iteration, exploration trees are supplied as a tuple with source to target as the first"""

        self.source_target = ETree.create(source=settings.source, target=settings.target)
        """:class:`~molpher.core.ExplorationTree.ExplorationTree` from source to target"""
        self.target_source = ETree.create(source=settings.target, target=settings.source)
        """:class:`~molpher.core.ExplorationTree.ExplorationTree` from target to source"""

        if settings.tree_params:
            self.source_target.params = settings.tree_params
            self.target_source.params = settings.tree_params

        self.source_target.thread_count = settings.max_threads
        self.target_source.thread_count = settings.max_threads

        self._iteration = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper(settings.filters)
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

        self.path = []
        """a list of SMILES strings representing the found path (defaults to an empty `list`)"""


    def find_minima(self):
        source_target_min = FindClosest()
        target_source_min = FindClosest()
        if self.verbose:
            print('Traversal times:')

            source_target_time = timeit(lambda : self.source_target.traverse(source_target_min))
            print('\tsource -> target: {0}'.format(source_target_time))
            target_source_time = timeit(lambda : self.target_source.traverse(target_source_min))
            print('\ttarget -> source: {0}'.format(target_source_time))

            print('\ttotal execution time: {0}'.format(source_target_time + target_source_time))

            print('Current Targets:')
            print('\tsource to target:', self.source_target.params['target'])
            print('\ttarget to source:', self.target_source.params['target'])
        else:
            self.source_target.traverse(source_target_min)
            self.target_source.traverse(target_source_min)

        return source_target_min.closest, target_source_min.closest

    def __call__(self):
        """
        Execute the search and return the path.

        :return: a list of SMILES strings representing the found path (defaults to an empty `list`)
        """

        counter = 0
        connecting_molecule = None
        while True:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self._iteration:
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

            if self.iter_callback:
                self.iter_callback((self.source_target, self.target_source,))

            source_target_min, target_source_min = self.find_minima()

            print('Current Minima:')
            print('\tsource to target:', source_target_min.getSMILES(), source_target_min.getDistToTarget())
            print('\ttarget to source:', target_source_min.getSMILES(), target_source_min.getDistToTarget())

            self.source_target.params = {
                'target' : target_source_min.getSMILES()
            }
            self.target_source.params = {
                'target' : source_target_min.getSMILES()
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

        source_target_path = find_path(self.source_target, connecting_molecule)
        target_source_path = find_path(self.target_source, connecting_molecule)
        assert source_target_path.pop(-1) == connecting_molecule
        target_source_path.reverse()
        source_target_path.extend(target_source_path)
        self.path = source_target_path
        print('Path:', self.path)
        return self.path