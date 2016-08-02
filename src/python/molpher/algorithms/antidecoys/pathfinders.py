from statistics import mean

from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import ExtendTreeOper
from molpher.core.operations import FilterMorphsOper
from molpher.core.operations import GenerateMorphsOper
from molpher.core.operations import PruneTreeOper
from molpher.core.operations import SortMorphsOper
from molpher.core.operations import CleanMorphsOper
from molpher.core.operations.callbacks import SortMorphsCallback

from .utils import timeit, evaluate_path
from .custom_opers import TopXFilter, GatherAntiFPScores

from .settings import MAX_THREADS, MAX_ITERS_PER_PATH, WAIT_FOR_ANTIDECOYS, ANTIDECOYS_DISTANCE_SWITCH, ROLLBACK_MAX_ITERS, RESET_CLOSEST_THRESHOLD, ROLLBACK_MAX_ITERS_ON_CLOSEST, MAX_ANTIFP_SURVIVORS, COMMON_BITS_PERC_THRS, MIN_ANTIDECOY_ITERS

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

    class FindTopClosest:

        def __init__(self, threshold):
            self.top_closest = []
            self.threshold = threshold

        def __call__(self, morph):
            if morph.dist_to_target < self.threshold:
                self.top_closest.append(morph.smiles)

    class AntiFpSortCallback(SortMorphsCallback):

        def __init__(self, antifp_scores):
            super(BidirectionalPathFinder.AntiFpSortCallback, self).__init__()
            self.minimum_common_bits_perc = 1.0
            self.antifp_scores = antifp_scores

        def __call__(self, a, b):
            perc_a = self.antifp_scores[a.getSMILES()]
            perc_b = self.antifp_scores[b.getSMILES()]
            minimum = min(perc_a, perc_b, self.minimum_common_bits_perc)
            if minimum < self.minimum_common_bits_perc:
                self.minimum_common_bits_perc = minimum
            return perc_a < perc_b

    def __init__(self, source, target, verbose=True, antidecoys_filter=None, paths_antifingerprint=None, antifingerprint=None):
        self.source = source
        self.target = target
        self.antidecoys_filter = antidecoys_filter
        self.paths_antifingerprint = paths_antifingerprint
        self.verbose = verbose

        self.source_target = ETree.create(source=source, target=target)
        self.source_target.thread_count = MAX_THREADS

        self.target_source = ETree.create(source=target, target=source)
        self.target_source.thread_count = MAX_THREADS

        self.source_target_min = self.FindClosest()
        self.target_source_min = self.FindClosest()

        if self.verbose:
            print("Tree Parameters:")
            print('\tsource -> target: {0}'.format(self.source_target.params))
            print('\ttarget -> source: {0}'.format(self.target_source.params))

        self._antifp_scores = None
        self._antifp_sort_callback = None
        self.antifingerprint = None
        if antifingerprint:
            self._antifp_scores = dict()
            self._antifp_sort_callback = self.AntiFpSortCallback(antifp_scores=self._antifp_scores)
            self.antifingerprint = antifingerprint

        if self._antifp_sort_callback:
            self._iteration = [
                GenerateMorphsOper()
                , FilterMorphsOper(FilterMorphsOper.SYNTHESIS | FilterMorphsOper.WEIGHT | FilterMorphsOper.DUPLICATES | FilterMorphsOper.HISTORIC_DESCENDENTS | FilterMorphsOper.MAX_DERIVATIONS, self.verbose)
                , CleanMorphsOper() if self._antifp_sort_callback else None
                , GatherAntiFPScores(self._antifp_scores, self.antifingerprint)
                , SortMorphsOper(callback=self._antifp_sort_callback) if self._antifp_sort_callback else None
                , TopXFilter(MAX_ANTIFP_SURVIVORS) if self._antifp_sort_callback else None
                , CleanMorphsOper()
                , SortMorphsOper()
                , FilterMorphsOper(FilterMorphsOper.PROBABILITY, self.verbose)
                , self.antidecoys_filter
                , ExtendTreeOper()
                , PruneTreeOper()
            ]
        else:
            self._iteration = [
                GenerateMorphsOper()
                , SortMorphsOper()
                , FilterMorphsOper(self.verbose)
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

    @staticmethod
    def fetch_top_closest(tree, threshold=RESET_CLOSEST_THRESHOLD):
        find_top_closest = BidirectionalPathFinder.FindTopClosest(threshold)
        tree.traverse(find_top_closest)
        return find_top_closest.top_closest
        # return [x.smiles for x in tree.leaves if x.dist_to_target < threshold]

    def _rollback_closest(self, tree, rollback_max_iters):
        for x in self.fetch_top_closest(tree):
            if tree.hasMol(x):
                self._rollback_path(tree, start_mol=x, rollback_max_iters=rollback_max_iters)

    def reset(self, connecting_molecule=None, max_iters_rollback=ROLLBACK_MAX_ITERS, max_iters_rollback_closest=ROLLBACK_MAX_ITERS_ON_CLOSEST):
        if not connecting_molecule:
            assert self.connecting_molecule
            connecting_molecule = self.connecting_molecule
        self._rollback_path(tree=self.source_target, start_mol=connecting_molecule, rollback_max_iters=max_iters_rollback)
        self._rollback_path(tree=self.target_source, start_mol=connecting_molecule, rollback_max_iters=max_iters_rollback)

        self._rollback_closest(tree=self.source_target, rollback_max_iters=max_iters_rollback_closest)
        self._rollback_closest(tree=self.target_source, rollback_max_iters=max_iters_rollback_closest)

        self.update_target(tree=self.source_target, target=self.target)
        self.update_target(tree=self.target_source, target=self.source)
        self.source_target_min = self.FindClosest()
        self.target_source_min = self.FindClosest()

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
                if oper.__class__ == TopXFilter:
                    print("Top 30 morphs (antidecoys):")
                    source_target_mask = self.source_target.candidates_mask
                    target_source_mask = self.target_source.candidates_mask
                    source_target_candidates = self.source_target.candidates
                    target_source_candidates = self.target_source.candidates
                    source_target_mins = [self._antifp_scores[mol.smiles] for idx, mol in enumerate(source_target_candidates) if source_target_mask[idx]]
                    target_source_mins = [self._antifp_scores[mol.smiles] for idx, mol in enumerate(target_source_candidates) if target_source_mask[idx]]
                    print('\tsource -> target: {0}'.format(source_target_mins[:30]))
                    print('\ttarget -> source: {0}'.format(target_source_mins[:30]))
                    scores = source_target_mins + target_source_mins
                    mean_score = mean(scores)
                    print("Mean antidecoys score: {0}".format(mean_score))
                    if mean_score < COMMON_BITS_PERC_THRS:
                        antidecoys_off = True
                        print("Antidecoys turned off (COMMON_BITS_PERC_THRS).")
                    print("Antifingerprint distance (absolute minimum):", self._antifp_sort_callback.minimum_common_bits_perc)

            if antidecoys_off and counter >= MIN_ANTIDECOY_ITERS:
                self._iteration = [
                    GenerateMorphsOper()
                    , SortMorphsOper()
                    , FilterMorphsOper(self.verbose)
                    , ExtendTreeOper()
                    , PruneTreeOper()
                ]
            if self._antifp_scores:
                self._antifp_scores.clear()

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
                print("Antidecoys turned off (ANTIDECOYS_DISTANCE_SWITCH).")

            self.update_target(self.source_target, self.target_source_min.closest.getSMILES())
            self.target_source_min = self.FindClosest()
            self.update_target(self.target_source, self.source_target_min.closest.getSMILES())
            self.source_target_min = self.FindClosest()

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
                if self.paths_antifingerprint:
                    path_valid, common_perc = evaluate_path(self.path, self.paths_antifingerprint)
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
            self.path = None
            return None