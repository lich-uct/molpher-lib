from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import ExtendTreeOper
from molpher.core.operations import FilterMorphsOper
from molpher.core.operations import GenerateMorphsOper
from molpher.core.operations import PruneTreeOper
from molpher.core.operations import SortMorphsOper

from .utils import timeit, evaluate_path

from .settings import MAX_THREADS, MAX_ITERS_PER_PATH, WAIT_FOR_ANTIDECOYS, ANTIDECOYS_DISTANCE_SWITCH, ROLLBACK_MAX_ITERS, RESET_CLOSE_THRESHOLD

class BidirectionalPathFinder:

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

    def __init__(self, source, target, verbose=True, antidecoys_filter=None, path_antifingerprint=None):
        self.source = source
        self.target = target
        self.antidecoys_filter = antidecoys_filter
        self.path_antifingerprint = path_antifingerprint
        self.verbose = verbose

        self.source_target = ETree.create(source=source, target=target)
        self.source_target.thread_count = MAX_THREADS

        self.target_source = ETree.create(source=target, target=source)
        self.target_source.thread_count = MAX_THREADS

        self.source_target_min = self.FindClosest()
        self.target_source_min = self.FindClosest()

        print("Tree Parameters:")
        print('\tsource -> target: {0}'.format(self.source_target.params))
        print('\ttarget -> source: {0}'.format(self.target_source.params))

        self._iteration = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper(self.verbose)
            , self.antidecoys_filter
            , ExtendTreeOper()
            , PruneTreeOper()
        ]
        self._iteration = [x for x in self._iteration if x]
        self.path = []
        self.connecting_molecule = None

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

    @staticmethod
    def find_top_parent(tree, smile, iters):
        counter = 0
        next_smile = smile
        previous_smile = smile
        while counter < iters:
            current = tree.fetchMol(next_smile)
            next_smile = current.getParentSMILES()
            if not next_smile:
                return previous_smile
            previous_smile = current.getSMILES()
            counter += 1
        return previous_smile

    def _rollback_path(self, tree, start_mol, rollback_max_iters):
        top_smile = self.find_top_parent(tree, start_mol, rollback_max_iters)
        tree.deleteSubtree(top_smile)
        assert not tree.hasMol(start_mol)

    def update_target(self, tree, target):
        if target != tree.params['source']:
            tree.params = {
                'target' : target
            }
        self.source_target_min = self.FindClosest()
        self.target_source_min = self.FindClosest()


    @staticmethod
    def fetch_top_closest(tree, threshold=RESET_CLOSE_THRESHOLD):
        return [x.smiles for x in tree.leaves if x.dist_to_target < threshold]

    def _rollback_closest(self, tree, rollback_max_iters):
        for x in self.fetch_top_closest(tree):
            if tree.hasMol(x):
                self._rollback_path(tree, start_mol=x, rollback_max_iters=rollback_max_iters)

    def reset(self, connecting_molecule=None, max_iters_rollback=ROLLBACK_MAX_ITERS):
        if not connecting_molecule:
            assert self.connecting_molecule
            connecting_molecule = self.connecting_molecule
        self._rollback_path(tree=self.source_target, start_mol=connecting_molecule, rollback_max_iters=max_iters_rollback)
        self._rollback_path(tree=self.target_source, start_mol=connecting_molecule, rollback_max_iters=max_iters_rollback)

        self._rollback_closest(tree=self.source_target, rollback_max_iters=max_iters_rollback)
        self._rollback_closest(tree=self.target_source, rollback_max_iters=max_iters_rollback)

        self.update_target(tree=self.source_target, target=self.target)
        self.update_target(tree=self.target_source, target=self.source)

        self.path = []
        self.connecting_molecule = None

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
                    print("Generated morphs:")
                    print('\tsource -> target: {0}'.format(len(self.source_target.candidates)))
                    print('\ttarget -> source: {0}'.format(len(self.target_source.candidates)))

            print("Accepted morphs:")
            print('\tsource -> target: {0}'.format(len(self.source_target.leaves)))
            print('\ttarget -> source: {0}'.format(len(self.target_source.leaves)))

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
            if min(source_target_min_dist, target_source_min_dist) < ANTIDECOYS_DISTANCE_SWITCH and self.antidecoys_filter:
                antidecoys_off = True
                print("Antidecoys turned off.")

            self.update_target(self.source_target, self.target_source_min.closest.getSMILES())
            self.update_target(self.target_source, self.source_target_min.closest.getSMILES())

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
            if self.target_source.path_found:
                connecting_molecule = self.target_source.params['target']
                if self.verbose:
                    print('Path Found in tree going from target to source')
                    print('Connecting molecule:', connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
            if connecting_molecule:
                source_target_path = self._find_path(self.source_target, connecting_molecule)
                target_source_path = self._find_path(self.target_source, connecting_molecule)
                assert source_target_path.pop(-1) == connecting_molecule
                target_source_path.reverse()
                source_target_path.extend(target_source_path)
                self.path = source_target_path

            if self.path:
                path_valid = True
                common_perc = None
                if self.path_antifingerprint:
                    path_valid, common_perc = evaluate_path(self.path, self.path_antifingerprint)
                if not path_valid:
                    print('Path will be removed due to antidecoys (% in common {0}): {1}'.format(common_perc, self.path))
                    print('Rolling back...')

                    self.reset(connecting_molecule)
                    connecting_molecule = None

                    continue
                else:
                    break

        if not search_failed:
            print('Path found:', self.path)
            self.connecting_molecule = connecting_molecule
            return self.path
        else:
            print('Search reached maximum number of iterations. Aborting...')
            return None