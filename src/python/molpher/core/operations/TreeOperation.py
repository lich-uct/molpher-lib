import molpher
import molpher.core.ExplorationTree

from abc import ABCMeta, abstractmethod

from molpher.core._utils import shorten_repr


class TreeOperation(molpher.swig_wrappers.core.TreeOperation):
    """
    Abstract base class derived from the `molpher.swig_wrappers.core.TreeOperation`
    proxy class.

    """

    __metaclass__ = ABCMeta

    def __repr__(self):
        return shorten_repr(TreeOperation, self)

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
            tree.__class__ = molpher.core.ExplorationTree # 'cast' the wrapped class to the 'pretty' Python proxy class
        return tree

    @abstractmethod
    def __call__(self):
        pass