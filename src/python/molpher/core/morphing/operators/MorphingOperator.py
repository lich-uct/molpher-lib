
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
from abc import ABCMeta, abstractmethod
from molpher.core._utils import shorten_repr


class MorphingOperator(molpher.swig_wrappers.core.MorphingOperator):
    """
    Abstract base class derived from the `molpher.swig_wrappers.core.MorphingOperator`
    proxy class.

    """

    __metaclass__ = ABCMeta

    def __init__(self):
        super(MorphingOperator, self).__init__()
        self.original = self.getOriginal()

    def __repr__(self):
        return shorten_repr(MorphingOperator, self)

    @abstractmethod
    def setOriginal(self, mol):
        pass

    @abstractmethod
    def morph(self):
        pass