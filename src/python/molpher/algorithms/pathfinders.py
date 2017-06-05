from molpher.algorithms.functions import find_path
from molpher.core import ExplorationTree as ETree


class BasicPathfinder:
    """
    :param settings: settings to use in the search
    :type settings: `Settings`

    A very basic pathfinder class that can be used to run exploration with
    any combination of operations.
    """

    def __init__(self, settings, operations):
        self.settings = settings
        """a settings class (should be a subclass of `Settings`)"""

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
            if counter > self.settings.max_iters:
                break
            counter+=1
            print('Iteration {0}'.format(counter))
            for oper in self._iteration:
                self.tree.runOperation(oper)

        self.path = find_path(self.tree, self.tree.params['target'])
        print('Path found:', self.path)
        return self.path