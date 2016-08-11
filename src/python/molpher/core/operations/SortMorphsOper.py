
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

    class SortCallback(SortMorphsCallback):

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