
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