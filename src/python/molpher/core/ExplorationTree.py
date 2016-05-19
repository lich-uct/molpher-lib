
# Copyright (c) 2016 Martin Sicho
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import molpher
import warnings

from molpher.core.ExplorationData import ExplorationData
from molpher.core.MolpherMol import MolpherMol
from molpher.core._utils import shorten_repr
from molpher.core.operations import TraverseOper
from molpher.core.operations.callbacks import TraverseCallback

class Callback(TraverseCallback):
    """
    :param callback: the callable to call every time a molecule is encountered during traversal
    :type callback: any callable object with one required parameter

    Basic callback class used to traverse the tree with the `ExplorationTree.traverse()` method.

    It registers a callable and calls it every time a :term:`morph` is processed.

    """

    def __init__(self, callback):
        super(Callback, self).__init__()
        self._callback = callback

    def __call__(self, morph):
        """
        This method is called by the C++ code. It
        just calls the registered callback for a :term:`morph`

        :param morph: morph in the currently processed tree
        :type morph: `molpher.swig_wrappers.core.MolpherMol`
        """

        morph.__class__ = MolpherMol
        self._callback(morph.tree.fetchMol(morph.smiles)) # needs to be done this way (SEGFAULT otherwise, the input parameter seems to be destroyed by SWIG after the call)

class ExplorationTree(molpher.swig_wrappers.core.ExplorationTree):
    """
    This a specialized version of the `molpher.swig_wrappers.core.ExplorationTree` proxy class.
    It implements some additional functionality for ease of use from Python.

    ..  attention:: This class has no constructor defined. Use the :meth:`create`
            factory method to obtain instances of this class.

    .. seealso:: `molpher.swig_wrappers.core.ExplorationTree`

    """

    def __repr__(self):
        return shorten_repr(ExplorationTree, self)

    def __init__(self):
        super(ExplorationTree, self).__init__()

    @staticmethod
    def _cast_mols(mols):
        ret = [ x for x in mols ]
        for morph in ret:
            morph.__class__ = MolpherMol
        return ret

    @staticmethod
    def create(tree_data=None, source=None, target=None, callback_class=Callback):
        """
        Create an exploration tree.

        :param tree_data: the morphing parameters (optional if ``source`` and ``target`` are specified)
        :type tree_data: `molpher.swig_wrappers.core.ExplorationData` (or its derived class),
            a `dict` of parameters (in the same format as in :class:`~molpher.core.ExplorationData.ExplorationData`'s constructor)
            or a path to a :term:`XML template` or a :term:`tree snapshot`
        :param source: SMILES of the source molecule or a molecule directly (uses a copy not the instance itself)
        :type source: `str` or :class:`~molpher.core.MolpherMol.MolpherMol`
        :param target: SMILES of the target molecule
        :type target: `str`
        :param callback_class: the class to use when making a callback using the :meth:`traverse` method
        :type callback_class: any class derived from `Callback`

        ..  note:: When ``tree_data`` is specified, ``source`` and ``target`` are always ignored.

        """

        ret = None
        if tree_data and source and target:
            warnings.warn(
                "Both tree_data and source and target specified. Using the values in parameters..."
                , RuntimeWarning
            )
        if tree_data and (isinstance(tree_data, molpher.wrappers.ExplorationData) or type(tree_data) == str):
            ret = super(ExplorationTree, ExplorationTree).create(tree_data)
        elif tree_data:
            _params = ExplorationData(**tree_data)
            ret = super(ExplorationTree, ExplorationTree).create(_params)
        elif source and target:
            if type(source) == str or type(target) == str:
                if type(source) == MolpherMol:
                    source = source.getSMILES()
                if type(target) == MolpherMol:
                    target = target.getSMILES()
            ret = super(ExplorationTree, ExplorationTree).create(source, target)
        else:
            raise AttributeError('Invalid set of parameters specified.')

        if not ret:
            raise RuntimeError('No tree initilized.')

        ret.callback_class = callback_class
        ret.__class__ = ExplorationTree

        return ret

    @property
    def params(self):
        """
        A dictionary representing the current :term:`exploration parameters`.

        It is possible to assign a new dictionary (or an instance of
        the `molpher.swig_wrappers.core.ExplorationData` class)
        to update the current parameters.

        .. note:: Only parameters defined in the supplied dictionary are changed
            and if an instance of `molpher.swig_wrappers.core.ExplorationData`
            is supplied only the parameters are read from it (the tree structure
            remains the same).

        :return: current parameters
        :rtype: `dict`
        """
        data = ExplorationData(other=self.asData())
        return data.param_dict

    @params.setter
    def params(self, params):
        if isinstance(params, molpher.wrappers.ExplorationData):
            self.update(params)
        else:
            new_params = self.params
            new_params.update(params)
            data = ExplorationData(**new_params)
            self.update(data)

    @property
    def generation_count(self):
        """
        :return: Number of :term:`morph generations <morph generation>` connected to the tree so far.
        :rtype: `int`
        """

        return self.getGenerationCount()

    @property
    def path_found(self):
        """
        :return: `True` if the :term:`target molecule` is present in the tree, `False` otherwise.
        :rtype: `bool`
        """

        return self.isPathFound()

    @property
    def leaves(self):
        """
        :return: the current leaves of the tree
        :rtype: `tuple` of :class:`~molpher.core.MolpherMol.MolpherMol` instances
        """

        return tuple(self._cast_mols(self.fetchLeaves()))

    @property
    def candidates(self):
        """
        :return: the :term:`candidate morphs` (:term:`morphs <morph>` generated by a single call to `generateMorphs()`.)
        :rtype: `tuple` of :class:`~molpher.core.MolpherMol.MolpherMol` instances
        """
        return tuple(self._cast_mols(self.getCandidateMorphs()))

    @property
    def candidates_mask(self):
        """
        A `tuple` of `bool` objects that serve as means of filtering the :term:`candidate morphs`.
        Each :term:`morph` in `candidates` has a `bool` variable assigned to it in this `tuple`
        -- only morphs with `True` at the appropriate position are added to the tree when :meth:`extend`
        is called.

        It can be changed by assigning a new `tuple` or a call to `setCandidateMorphsMask()`.

        :return: currently selected :term:`candidate morphs` represented as a `tuple` of `bool` objects
        :rtype: `tuple`
        """

        return self.getCandidateMorphsMask()

    @candidates_mask.setter
    def candidates_mask(self, mask):
        self.setCandidateMorphsMask(mask)

    @property
    def thread_count(self):
        """
        :return: maximum number of threads this instance will use
        :rtype: `int`
        """

        return self.getThreadCount()

    @thread_count.setter
    def thread_count(self, val):
        self.setThreadCount(val)

    def fetchMol(self, canonSMILES):
        """
        Returns a molecule from the tree using a canonical
        SMILES string.

        Raises a `RuntimeError` if the molecule is not found.

        :param canonSMILES: SMILES string of the molecule to fetch
        :type canonSMILES: `str`
        :return: molecule from a tree
        :rtype: :class:`~molpher.core.MolpherMol.MolpherMol`
        """

        ret = super(ExplorationTree, self).fetchMol(canonSMILES)
        ret.__class__ = MolpherMol
        return ret

    def traverse(self, callback, start_mol = None):
        """
        This method can be used to traverse the whole tree structure (or just a subtree)
        starting from the root to leaves. It takes a callback function that accepts a single required argument
        and traverses the whole tree starting from its root (or root of a specified subtree -- see ``start_mol``) and
        calls the supplied callback with with encountered morph as its parameter.

        :param callback: the callback to call
        :type callback: a callable object that takes a single argument
        :param start_mol: the root of a subtree to explore as canonical SMILES
            or :py:class:`~molpher.core.MolpherMol.MolpherMol` instance
        :type start_mol: `str` or :py:class:`~molpher.core.MolpherMol.MolpherMol`
        """

        if start_mol and type(start_mol) == MolpherMol:
            TraverseOper(self, self.callback_class(callback), start_mol)()
        elif start_mol:
            mol = self.fetchMol(start_mol)
            TraverseOper(self, self.callback_class(callback), mol)()
        else:
            cb = self.callback_class(callback)
            TraverseOper(self, cb)()

    def asData(self):
        """
        :return: the tree as an :class:`~molpher.core.ExplorationData.ExplorationData` instance
        :rtype: :class:`~molpher.core.ExplorationData.ExplorationData`
        """

        return ExplorationData(other=super(ExplorationTree, self).asData())