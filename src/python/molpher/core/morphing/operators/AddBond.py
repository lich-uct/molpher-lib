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

import molpher
from molpher.core.MolpherMol import MolpherMol
from molpher.core._utils import shorten_repr

class AddBond(molpher.swig_wrappers.core.AddBond):

    def __repr__(self):
        return shorten_repr(self.__class__, self)

    def getOpenBonds(self):
        return super(AddBond, self).getOpenBonds()

    def morph(self):
        ret = super(AddBond, self).morph()
        if ret:
            ret.__class__ = MolpherMol
        return ret

    def getOriginal(self):
        ret = super(AddBond, self).getOriginal()
        if ret:
            ret.__class__ = MolpherMol
        return ret