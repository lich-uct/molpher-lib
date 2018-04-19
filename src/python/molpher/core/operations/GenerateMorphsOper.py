
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
from molpher.core.morphing.MorphCollector import MorphCollector


class GenerateMorphsOper(core.GenerateMorphsOper, TreeOperation):

    def __init__(self, tree = None, collectors=()):
        self._collectors = [MorphCollector(x) for x in collectors]
        if not self._collectors and not tree:
            super(GenerateMorphsOper, self).__init__()
        elif not self._collectors and tree:
            super(GenerateMorphsOper, self).__init__(tree)
        elif self._collectors and not tree:
            super(GenerateMorphsOper, self).__init__(self._collectors)
        elif self._collectors and tree:
            super(GenerateMorphsOper, self).__init__(tree, self._collectors)