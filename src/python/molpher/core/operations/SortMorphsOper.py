
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

import molpher.swig_wrappers.core
from molpher.core.operations.TreeOperation import TreeOperation
from molpher.core.operations.callbacks import SortMorphsCallback
from molpher.core.MolpherMol import MolpherMol


class SortMorphsOper(molpher.swig_wrappers.core.SortMorphsOper, TreeOperation):
    """
    :param tree: a tree to run this operation on (optional)
    :type tree: :class:`~molpher.core.ExplorationTree.ExplorationTree` instance
    :param callback: callable to call every time two molecules need to be compared (optional)
    :type callback: any callable object with two positional arguments and a single `bool` return value

    This operation sorts the `candidates` list of an :class:`~molpher.core.ExplorationTree.ExplorationTree` using a given callback function.
    If no callback function is specified, candidates are sorted in ascending order according to `dist_to_target`.

    """

    class SortCallback(SortMorphsCallback):
        """
        :param callback: callable to call every time two molecules need to be compared
        :type callback: any callable object with two positional arguments and a single `bool` return value

        Customized callback that calls any supplied
        callable as a callback of `SortMorphsOper`
        and injects instances of `molpher.core.MolpherMol`
        rather than :class:`molpher.swig_wrappers.core.MolpherMol`.

        """

        def __init__(self, callback):
            super(SortMorphsOper.SortCallback, self).__init__()
            self._callback = callback

        def __call__(self, a, b):
            a.__class__ = MolpherMol
            b.__class__ = MolpherMol
            return self._callback(a, b)


    def __init__(self, tree=None, callback=None):
        self._callback = None
        if not callback:
            self._callback = molpher.swig_wrappers.core.DefaultSortCallback()
        else:
            self._callback = self.SortCallback(callback)

        if tree:
            super(SortMorphsOper, self).__init__(tree, self._callback)
        else:
            super(SortMorphsOper, self).__init__(self._callback)