"""
General utility functions for the algorithms.
"""

import time

def timeit(func):
    """
    Executes :samp:`func` and
    returns its runtime.

    :param func: function to time
    :type func: any callable
    :return: runtime of :samp:`func`
    :rtype: `float`
    """

    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds


def find_path(tree, end_mol):
    """
    Backtracks the tree starting from :samp:`end_mol`
    and returns a list of SMILES strings representing
    the path from the root of the tree to :samp:`end_mol`.

    :param tree: a tree to backtrack through
    :type tree: instance of :class:`~molpher.core.ExplorationTree.ExplorationTree`
    :param end_mol: SMILES of the last molecule in the requested path
    :type end_mol: `str`
    :return: list of SMILES of molecules on the path
    :rtype: `list` of `str`
    """

    path = []
    current = tree.fetchMol(end_mol)
    path.append(current.getSMILES())
    while current != '':
        current = current.getParentSMILES()
        if current:
            current = tree.fetchMol(current)
            path.append(current.getSMILES())
    path.reverse()
    return path