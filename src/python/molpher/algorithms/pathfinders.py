from molpher.algorithms.functions import find_path
from molpher.core import ExplorationTree as ETree


class BasicPathfinder:
    """
    :param settings: settings to use in the search
    :type settings: `Settings`

    A very basic pathfinder class that can be used to run exploration with
    any combination of operations.
    """

    class MaxItersReachedException(Exception):

        def __init__(self, tree):
            super(BasicPathfinder.MaxItersReachedException, self).__init__(
                "Maximum number of  iterations reached while searching "
                "for a path\n\t source: {0}\n\t target: {1}".format(tree.source, tree.target))


    def __init__(self, settings, operations, iter_callback = None):
        self.settings = settings
        """a settings class (should be a subclass of `Settings`)"""

        self._iter_callback = iter_callback
        """a callable to run after finishing one iteration, it receives the exploration tree as its parameter"""

        self.tree = ETree.create(source=self.settings.source, target=self.settings.target)
        """:class:`~molpher.core.ExplorationTree.ExplorationTree` used in the search"""
        if self.settings.tree_params:
            self.tree.params = self.settings.tree_params
        self.tree.thread_count = self.settings.max_threads

        self._iteration = operations

        self.path = None
        """a list of SMILES strings if a path was found, `None` otherwise"""

    def __call__(self):
        """
        Executes the search

        :return: discovered path
        :rtype: `list` of `str`
        """

        counter = 0
        while not self.tree.path_found:
            counter+=1
            if counter > self.settings.max_iters:
                raise BasicPathfinder.MaxItersReachedException(self.tree)
            print('Iteration {0}'.format(counter))
            for oper in self._iteration:
                self.tree.runOperation(oper)

            if self._iter_callback:
                self._iter_callback(self.tree)

        self.path = find_path(self.tree, self.tree.params['target'])
        print('Path found:', self.path)
        return self.path