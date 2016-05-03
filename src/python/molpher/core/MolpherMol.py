import warnings

import molpher.swig_wrappers.core as wrappers
import molpher.core.ExplorationTree
from molpher.core._utils import shorten_repr


class MolpherMol(wrappers.MolpherMol):
    """
    :param smiles: smiles of the molecule that is to be created
    :type smiles: `str`
    :param other: another instance, the new instance will be its copy (tree ownership is not transferred onto the copy)
    :type other: `molpher.swig_wrappers.core.MolpherMol` or its derived class

    This a specialized version of the `molpher.swig_wrappers.core.MolpherMol` proxy class.
    It implements some additional functionality for ease of use from Python.

    .. seealso:: `molpher.swig_wrappers.core.MolpherMol`

    """

    def __init__(self, smiles=None, other=None):
        if smiles and other:
            warnings.warn(
                "Both SMILES string and another MolpherMol instance specified. Using another instance to initialize..."
                , RuntimeWarning
            )
        if other and isinstance(other, wrappers.MolpherMol):
            super(MolpherMol, self).__init__(other)
        elif smiles and type(smiles) == str:
            super(MolpherMol, self).__init__(smiles)
        elif not smiles and not other:
            super(MolpherMol, self).__init__()
        else:
            raise AttributeError('Invalid argumetns supplied to the constructor.')

    def __repr__(self):
        return shorten_repr(MolpherMol, self)

    def __lt__(self, other):
        return self.smiles < other.smiles

    def copy(self):
        """
        Returns a copy of this instance.
        If this instance has a `tree` assigned,
        the returned will have `None` instead.

        :return: a copy of this instance
        :rtype:
        """
        copy = super(MolpherMol, self).copy()
        copy.__class__ = MolpherMol
        return copy

    @property
    def tree(self):
        """
        A reference to the tree this instance is currently in.
        If the molecule is not present in any tree,
        this value is `None`.

        :return: reference to the tree this instance is currently in
        :rtype: :class:`~molpher.core.ExplorationTree.ExplorationTree` or `None`
        """

        ret = self.getTree()
        if ret:
            ret.__class__ = molpher.core.ExplorationTree
        return ret

    @property
    def smiles(self):
        """
        :return: canonical SMILES string of this molecule
        :rtype: `str`
        """

        return self.getSMILES()

    @smiles.setter
    def smiles(self, val):
        self.setSMILES(val)

    @property
    def parent_smiles(self):
        """
        Canonical SMILES string of the parent molecule
        in the tree.

        Can be an empty `str`, if the molecule is a
        root of the tree or is not associated with any.

        :return: canonical SMILES string of the parent molecule in the tree
        :rtype: `str`
        """

        return self.getParentSMILES()

    @property
    def parent_operator(self):
        """
        The name of the :term:`chemical operator <chemical operators>`
        :term:`selector <selectors>` that lead to the creation of this
        molecule.

        :return: name of the parent chemical operator
        :rtype: `str`
        """

        return wrappers.ChemOperShortDesc(self.getParentOper())

    @property
    def dist_to_target(self):
        """
        The value of the objective function.
        In the original implementation,
        this is the structural distance to the target
        molecule using a `similarity measure`.

        This value can be changed.

        :return: value of the objective function
        :rtype: `float`
        """

        return self.getDistToTarget()

    @dist_to_target.setter
    def dist_to_target(self, dist):
        self.setDistToTarget(dist)

    @property
    def sascore(self):
        """

        The synthetic feasibility score of the molecule according
        to Ertl.

        .. todo:: add reference

        This value can be changed.

        :return: synthetic feasibility score
        :rtype: `float`
        """

        return self.getSAScore()

    @sascore.setter
    def sascore(self, val):
        self.setSAScore(val)

    @property
    def historic_descendents(self):
        """
        Canonical SMILES strings of all molecules derived from
        this compound.

        :return:
        :rtype: `str`
        """
        return self.getHistoricDescendants()

    @property
    def descendents(self):
        """
        Canonical SMILES strings of all molecules derived from
        this compound that are currently present in the tree.

        :return:
        :rtype: `str`
        """
        return self.getDescendants()

    @property
    def gens_without_improvement(self):
        """
        Number of morph generations derived from this molecule that did not
        contain any morphs with an improvement in the objective function
        from the target molecule.

        This value can be changed.

        :return: number of non-producing generations
        :rtype: `int`
        """

        return self.getItersWithoutDistImprovement()

    @gens_without_improvement.setter
    def gens_without_improvement(self, value):
        self.setItersWithoutDistImprovement(value)