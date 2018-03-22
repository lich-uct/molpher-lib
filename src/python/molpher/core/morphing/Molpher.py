
# Copyright (c) 2017 Martin Sicho
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

import molpher.swig_wrappers.core as wrappers
from molpher.core.MolpherMol import MolpherMol
from molpher.core._utils import shorten_repr


class Molpher(wrappers.Molpher):
    """
    :param molecule: The starting molecule.
    :type molecule: an instance of 'molpher.swig_wrappers.core.MolpherMol' or its derived class.

    This a specialized version of the `molpher.swig_wrappers.core.morphing.Molpher` proxy class.
    It implements some additional functionality for ease of use from Python.

    .. seealso:: `molpher.swig_wrappers.core.AtomLibrary`
    """

    def __repr__(self):
        return shorten_repr(self.__class__, self)

    def __init__(self, molecule, operators, threads=0, attempts=30, max_iters=None):
        super(Molpher, self).__init__(molecule, operators, threads, attempts)
        self._iter_morphs_cache = []
        self._iter_counter = 0
        self.max_iters = max_iters

    def getMorphs(self):
        ret = super(Molpher, self).getMorphs()
        for x in ret:
            x.__class__ = MolpherMol
        return ret

    def __iter__(self):
        return self

    def __next__(self):
        if self.max_iters and self.max_iters > 0:
            self._iter_counter+=1
            if self._iter_counter > self.max_iters:
                raise StopIteration
        if self._iter_morphs_cache:
            return self._iter_morphs_cache.pop()
        else:
            self()
            self._iter_morphs_cache = list(self.getMorphs())
            return self._iter_morphs_cache.pop()

    def next(self):
        return self.__next__()