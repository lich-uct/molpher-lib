from molpher.algorithms.operations import FindClosest
from molpher.algorithms.pathfinders import BasicPathfinder
from molpher.core.operations import *

class ClassicPathFinder(BasicPathfinder):
    """
    :param settings: settings to use in the search
    :type settings: `Settings`

    A callable class which implements the original molpher algorithm as published in [1]_.

    .. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_
    """

    def __init__(self, settings):
        super(ClassicPathFinder, self).__init__(settings, [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper(settings.verbose)
            , ExtendTreeOper()
            , PruneTreeOper()
        ])

        self.find_closest = FindClosest()
        """instance of `FindClosest` that holds the molecule currently closest to target"""