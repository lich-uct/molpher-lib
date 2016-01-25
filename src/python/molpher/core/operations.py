"""
This module contains references to
the built-in operations defined
in `molpher.swig_wrappers.core`
and defines a more convenient
`TreeOperation` base class.

"""

import molpher
from molpher.core import ExplorationTree
from molpher.swig_wrappers.core import \
    GenerateMorphsOper \
    , FilterMorphsOper \
    , ExtendTreeOper \
    , PruneTreeOper \
    , SortMorphsOper \
    , FindLeavesOper \
    , TraverseOper \
    , TraverseCallback

from abc import ABCMeta, abstractmethod


class TreeOperation(molpher.swig_wrappers.core.TreeOperation):
    __metaclass__ = ABCMeta

    @property
    def tree(self):
        """
        The :class:`~molpher.core.ExplorationTree` instance this operation
        currently operates on. It can also be written into to change the
        current tree instance.

        :return: current :class:`~molpher.core.ExplorationTree` instance
        :rtype: :class:`~molpher.core.ExplorationTree`
        """

        return self.getTree()

    @tree.setter
    def tree(self, tree):
        self.setTree(tree)

    def getTree(self):
        """
        Getter which returns the :class:`~molpher.core.ExplorationTree` instance this operation
        currently operates on.

        :return: current :class:`~molpher.core.ExplorationTree` instance
        :rtype: :class:`~molpher.core.ExplorationTree`
        """

        tree = super(TreeOperation, self).getTree()
        if tree:
            tree.__class__ = ExplorationTree # 'cast' the wrapped class to the 'pretty' Python proxy class
        return tree

    @abstractmethod
    def __call__(self):
        pass