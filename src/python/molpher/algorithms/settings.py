import os

import multiprocessing


class Settings:
    """
    Holds basic settings and parameters shared among the exploration algorithms.
    """

    def __init__(
            self
            , source
            , target
            , storage_dir
            , max_threads = None
            , max_iters = 100
            , tree_params = None
            , verbose = False
    ):
        self.verbose = verbose
        """some algorithms can support verbose output"""

        self.source = source
        """SMILES string of the source molecule"""
        self.target = target
        """SMILES string of the target molecule"""

        self.storage_dir = storage_dir
        """path to a directory where results will be stored"""
        if not os.path.exists(self.storage_dir):
            os.mkdir(self.storage_dir)

        self.max_threads = max_threads if max_threads not in (None, False) else multiprocessing.cpu_count()
        """maximum number of threads to use in computations (uses all available CPUs by default)"""
        self.max_iters = max_iters
        """maximum number of iterations to spend on a search"""
        self.tree_params = tree_params
        """parameters of the exploration tree; either a dictionary like in `params` or an instance of :py:class:`~molpher.core.ExplorationData.ExplorationData`"""