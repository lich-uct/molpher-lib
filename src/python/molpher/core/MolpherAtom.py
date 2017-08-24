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
from molpher.core._utils import shorten_repr

class MolpherAtom(wrappers.MolpherAtom):
    """
    :param symbol: The symbol in the periodic table of the atom to create.
    :type symbol: `str`
    :param formal_charge: its formal charge (zero by default)
    :type formal_charge: `int`

    This a specialized version of the `molpher.swig_wrappers.core.MolpherAtom` proxy class.
    It implements some additional functionality for ease of use from Python.

    .. seealso:: `molpher.swig_wrappers.core.MolpherAtom`

    """

    def __init__(self, symbol, formal_charge=0):
        super(MolpherAtom, self).__init__(symbol, formal_charge)

    def __repr__(self):
        return shorten_repr(MolpherAtom, self)

    @property
    def symbol(self):
        return self.getSymbol()

    @property
    def mass(self):
        return self.getMass()

    @property
    def atomic_number(self):
        return self.getAtomicNum()

    @property
    def formal_charge(self):
        return self.getFormalCharge()

    @formal_charge.setter
    def formal_charge(self, val):
        self.setFormalCharge(val)

    @property
    def locking_mask(self):
        return self.getLockingMask()

    @locking_mask.setter
    def locking_mask(self, mask):
        self.setLockingMask(mask)

    @property
    def is_locked(self):
        return self.isLocked()

    @property
    def lock_info(self):
        return {
            'UNLOCKED' : not self.is_locked
            , 'NO_MUTATION' : bool(self.locking_mask & MolpherAtom.NO_MUTATION)
            , 'NO_ADDITION' : bool(self.locking_mask & MolpherAtom.NO_ADDITION)
            , 'KEEP_NEIGHBORS' : bool(self.locking_mask & MolpherAtom.KEEP_NEIGHBORS)
            , 'FULL_LOCK' : bool((self.locking_mask & MolpherAtom.FULL_LOCK) == MolpherAtom.FULL_LOCK)
        }
