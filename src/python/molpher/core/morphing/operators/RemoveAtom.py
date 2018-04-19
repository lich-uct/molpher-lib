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
from molpher.core.MolpherAtom import MolpherAtom
from molpher.core._utils import shorten_repr

class RemoveAtom(wrappers.RemoveAtom):
    """
    This a specialized version of the `molpher.swig_wrappers.core.morphing.operators.RemoveAtom` proxy class.
    It implements some additional functionality for ease of use from Python.

    .. seealso:: `molpher.swig_wrappers.core.RemoveAtom`

    """

    def __repr__(self):
        return shorten_repr(self.__class__, self)

    def getMarkedAtoms(self):
        ret = super(RemoveAtom, self).getMarkedAtoms()
        for x in ret:
            x.__class__ = MolpherAtom
        return ret

    def morph(self):
        ret = super(RemoveAtom, self).morph()
        if ret:
            ret.__class__ = MolpherMol
        return ret

    def getOriginal(self):
        ret = super(RemoveAtom, self).getOriginal()
        if ret:
            ret.__class__ = MolpherMol
        return ret