
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
from molpher.core.MolpherAtom import MolpherAtom
from molpher.core._utils import shorten_repr


class AtomLibrary(molpher.swig_wrappers.core.AtomLibrary):
    """
    :param atoms: A set of atoms that defines the library.
    :type atoms: an iterable of :class:`~.morphing.AtomLibrary.AtomLibrary` instances, another :class:`~.morphing.AtomLibrary.AtomLibrary` instance or simply a list of atom symbols.
    :type atom_probabilities: an iterable of numbers specifying the probability of each atom being incorporated in the generated structures. Needs to be the same length as :command:`atoms` or empty (equal probability for each atom).
    """

    def __init__(self, atoms, atom_probabilities=tuple()):
        atoms_ = []
        try:
            iter(atoms)
            if len(atoms) != 0:
                first = atoms[0]
                if type(first) == str:
                    for atom in atoms:
                        atoms_.append(MolpherAtom(atom))
            else:
                raise RuntimeError("Atom list must not be empty.")
        except TypeError:
            pass
        if (len(atoms) != len(atom_probabilities)) and not (len(atom_probabilities) == 0):
            raise RuntimeError("The length of the atom list does not match the length of the probability list.")
        super(AtomLibrary, self).__init__(atoms_ if atoms_ else atoms, atom_probabilities)

    def __repr__(self):
        return shorten_repr(self.__class__, self)

    @staticmethod
    def getDefaultLibrary():
        ret = molpher.swig_wrappers.core.AtomLibrary.getDefaultLibrary()
        ret.__class__ = AtomLibrary
        return ret

    def getAtoms(self):
        ret = super(AtomLibrary, self).getAtoms()
        for x in ret:
            x.__class__ = MolpherAtom
        return ret

    @property
    def atoms(self):
        return self.getAtoms()

    @property
    def atom_probabilities(self):
        return super(AtomLibrary, self).getAtomProbabilities()

    def getRandomAtom(self):
        ret = super(AtomLibrary, self).getRandomAtom()
        ret.__class__ = MolpherAtom
        return ret