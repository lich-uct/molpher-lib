
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

import molpher.swig_wrappers.core as core
from molpher.core.operations.TreeOperation import TreeOperation
from molpher.core.operations.callbacks import TraverseCallback
from molpher.core.MolpherMol import MolpherMol

class Callback(TraverseCallback):
    """
    :param callback: the callable to call every time a molecule is encountered during traversal
    :type callback: any callable object with one required parameter

    Basic callback class used to traverse the tree with an arbitrary callable.
    It registers a callable and calls it every time a :term:`morph` is processed.
    It also prevents possible errors and supplies instances of
    `molpher.core.MolpherMol` rather than `molpher.swig_wrappers.core.MolpherMol`.

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
        self._callback(morph.tree.fetchMol(
            morph.smiles))  # needs to be done this way (SEGFAULT otherwise, the input parameter seems to be destroyed by SWIG after the call)

class TraverseOper(core.TraverseOper, TreeOperation):

    def __init__(self, tree=None, callback=None, root=None):
        if not callback:
            raise RuntimeError('Callback must be specified.')
        self._callback = Callback(callback)

        if tree:
            super(TraverseOper, self).__init__(tree, self._callback)
        elif tree and root:
            super(TraverseOper, self).__init__(tree, self._callback, root)
        else:
            super(TraverseOper, self).__init__(self._callback)