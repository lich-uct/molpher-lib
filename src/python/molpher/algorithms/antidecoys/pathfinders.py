from statistics import mean

from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import ExtendTreeOper
from molpher.core.operations import FilterMorphsOper
from molpher.core.operations import GenerateMorphsOper
from molpher.core.operations import PruneTreeOper
from molpher.core.operations import SortMorphsOper
from molpher.core.operations import CleanMorphsOper
from molpher.core.operations.callbacks import SortMorphsCallback

from .utils import timeit, update_target, find_path
from .custom_opers import TopScoringFilter, GatherAntiFPScores

from .settings import AntidecoysSettings

class PathFinder:

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

    class AntiFpSortCallback(SortMorphsCallback):

        def __init__(self, antifp_scores):
            super(PathFinder.AntiFpSortCallback, self).__init__()
            self.minimum_common_bits_perc = 1.0
            self.maximum_common_bits_perc = 0.0
            self.antifp_scores = antifp_scores

        def __call__(self, a, b):
            perc_a = self.antifp_scores[a.getSMILES()]
            perc_b = self.antifp_scores[b.getSMILES()]
            minimum = min(perc_a, perc_b, self.minimum_common_bits_perc)
            maximum = max(perc_a, perc_b, self.minimum_common_bits_perc)
            if minimum < self.minimum_common_bits_perc:
                self.minimum_common_bits_perc = minimum
            if maximum > self.maximum_common_bits_perc:
                self.maximum_common_bits_perc = maximum
            return perc_a < perc_b

    def __init__(
            self
            , source
            , target
            , verbose=True
            , settings=AntidecoysSettings()
            , antifingerprint=None
    ):
        self.source = source
        self.target = target
        self.verbose = verbose
        self.settings = settings

        self.source_target = ETree.create(source=source, target=target)
        self.source_target.thread_count = self.settings.max_threads

        self.target_source = ETree.create(source=target, target=source)
        self.target_source.thread_count = self.settings.max_threads

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

        if self.antifingerprint:
            self._iteration = [
                GenerateMorphsOper()
                , FilterMorphsOper(
                    FilterMorphsOper.SYNTHESIS
                    | FilterMorphsOper.WEIGHT
                    | FilterMorphsOper.DUPLICATES
                    | FilterMorphsOper.HISTORIC_DESCENDENTS
                    | FilterMorphsOper.MAX_DERIVATIONS
                    , self.verbose
                )
                , CleanMorphsOper()
                , GatherAntiFPScores(
                    self._antifp_scores
                    , self.antifingerprint
                    , self.settings
                )
                , SortMorphsOper(callback=self._antifp_sort_callback)
                , TopScoringFilter(
                    self._antifp_scores
                    , self.settings.common_bits_max_thrs
                    , self.settings.min_accepted
                )
                , CleanMorphsOper()
                , SortMorphsOper()
                , FilterMorphsOper(FilterMorphsOper.PROBABILITY, self.verbose)
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

    def __call__(self):
        counter = 0
        connecting_molecule = None
        max_iters_reached = False
        antidecoys_off = False
        while True:
            counter+=1
            if counter > self.settings.max_iters_per_path:
                max_iters_reached = True
                break
            if not antidecoys_off and self.antifingerprint and counter > self.settings.max_iters:
                print("Maximum number of iterations with antidecoys reached ({0}).".format(self.settings.max_iters))
                antidecoys_off = True
            print('## Iteration {0} ##'.format(counter))
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
                if oper.__class__ == TopScoringFilter:
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
                    if not antidecoys_off and self.antifingerprint and mean_score < self.settings.common_bits_mean_thrs:
                        antidecoys_off = True
                        print("Mean antidecoys score threshold reached ({0}).".format(self.settings.common_bits_mean_thrs))

            if antidecoys_off and self.antifingerprint and counter >= self.settings.min_iters:
                print("Antidecoys turned off. Setting up algorithm for normal search...")
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
                print('Current Targets:')
                print('\tsource to target:', self.source_target.params['target'])
                print('\ttarget to source:', self.target_source.params['target'])

            self.source_target.traverse(self.source_target_min)
            self.target_source.traverse(self.target_source_min)

            source_target_min_dist = self.source_target_min.closest.getDistToTarget()
            target_source_min_dist = self.target_source_min.closest.getDistToTarget()
            print('Current Minima:')
            print('\tsource to target:', self.source_target_min.closest.getSMILES(), source_target_min_dist)
            print('\ttarget to source:', self.target_source_min.closest.getSMILES(), target_source_min_dist)
            if not antidecoys_off and min(source_target_min_dist, target_source_min_dist) < self.settings.distance_thrs:
                antidecoys_off = True
                print("Antidecoys turned off. Trees are sufficinetly close ({0}).".format(self.settings.distance_thrs))

            update_target(self.source_target, self.target_source_min.closest.getSMILES())
            # self.target_source_min = self.FindClosest()
            update_target(self.target_source, self.source_target_min.closest.getSMILES())
            # self.source_target_min = self.FindClosest()

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
                source_target_path = find_path(self.source_target, connecting_molecule)
                target_source_path = find_path(self.target_source, connecting_molecule)
                assert source_target_path.pop(-1) == connecting_molecule
                target_source_path.reverse()
                source_target_path.extend(target_source_path)
                self.path = source_target_path

            if self.path:
                break

        if not max_iters_reached:
            print('Path found:', self.path)
            self.connecting_molecule = connecting_molecule
            return self.path
        else:
            print('Search reached maximum number of iterations. Aborting...')
            self.path = None
            return None