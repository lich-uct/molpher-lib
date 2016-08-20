from molpher.algorithms.functions import find_path
from molpher.algorithms.operations import FindClosest
from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

class ClassicPathFinder:
    """
    :param settings: settings to use in the search
    :type settings: `Settings`

    A callable class which implements the original molpher algorithm as published in [1]_.

    .. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_
    """

    def __init__(self, settings):
        self.settings = settings
        """used settings as `Settings`"""

        self.tree = ETree.create(source=self.settings.source, target=self.settings.target)
        """:class:`~molpher.core.ExplorationTree.ExplorationTree` used in the search"""

        if self.settings.tree_params:
            self.tree.params = self.settings.tree_params
        self.tree.thread_count = self.settings.max_threads

        self._iteration = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper(settings.verbose)
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

        self.find_closest = FindClosest()
        """instance of `FindClosest` that holds the molecule currently closest to target"""

        self.path = None
        """if a path is found it is written here (defaults to `None`)"""

    def __call__(self):
        """
        Executes the search

        :return: found path
        :rtype: `list` of `str`
        """

        counter = 0
        while not self.tree.path_found:
            if counter > self.settings.max_iters:
                break
            counter+=1
            print('Iteration {0}'.format(counter))
            for oper in self._iteration:
                self.tree.runOperation(oper)
            self.tree.traverse(self.find_closest)
            print('Closest to target: {0} {1}'.format(self.find_closest.closest.smiles, self.find_closest.closest.dist_to_target))

        self.path = find_path(self.tree, self.tree.params['target'])
        print('Path found:', self.path)
        return self.path