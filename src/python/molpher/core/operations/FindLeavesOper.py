import molpher.swig_wrappers.core
from molpher.core.operations.TreeOperation import TreeOperation
from molpher.core.MolpherMol import MolpherMol

class FindLeavesOper(molpher.swig_wrappers.core.FindLeavesOper, TreeOperation):

    @property
    def leaves(self):
        """
        The last set of leaves fetched from a tree.

        :return: fetched leaves
        :rtype: `tuple` of :class:`~molpher.core.MolpherMol.MolpherMol` instances
        """

        ret = self.fetchLeaves()
        for mol in ret:
            mol.__class__ = MolpherMol
        return ret