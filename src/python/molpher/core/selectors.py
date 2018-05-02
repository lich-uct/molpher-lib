"""
Contains all global selectors that are usually
used when creating an :term:`exploration tree`
or setting any of its parameters during runtime.

There are three types of selectors:

    #. fingerprints selectors
            Their names are prepended with 'FP\_' and are used to
            either set the `fingerprint` member of the
            :class:`~molpher.core.ExplorationData.ExplorationData`
            class or as the value of the *fingerprint* key when calling
            :meth:`~molpher.core.ExplorationTree.ExplorationTree.create`
            with the ``params`` parameter or writing into the `params`
            member of the :class:`~molpher.core.ExplorationTree.ExplorationTree`.

            `FP_MORGAN` is the default option.

    #. similarity coefficient selectors
            Their names are prepended with 'SC\_' and are used to
            either set the `similarity` member of the
            :class:`~molpher.core.ExplorationData.ExplorationData`
            class or as the value of the *similarity* key when calling
            :meth:`~molpher.core.ExplorationTree.ExplorationTree.create`
            with the ``params`` parameter or writing into the `params`
            member of the :class:`~molpher.core.ExplorationTree.ExplorationTree`.

            `SC_TANIMOTO` is the default option.

    #. chemical operators
            Their names are prepended with 'OP\_' and an :term:`iterable`
            of them is used to
            either set the :attr:`~molpher.core.ExplorationData.ExplorationData.operators` member of the
            :class:`~molpher.core.ExplorationData.ExplorationData`
            class or as items of the :term:`iterable` assigned to *operators* key
            when calling :meth:`~molpher.core.ExplorationTree.ExplorationTree.create`
            with the ``params`` parameter or writing into the `params`
            member of the :class:`~molpher.core.ExplorationTree.ExplorationTree`.

            All of the available selectors are used by default.

"""


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

from molpher.swig_wrappers.core import \
    FP_ATOM_PAIRS\
    , FP_EXT_ATOM_PAIRS\
    , FP_EXT_MORGAN\
    , FP_EXT_TOPOLOGICAL\
    , FP_EXT_TOPOLOGICAL_LAYERED_1\
    , FP_EXT_TOPOLOGICAL_LAYERED_2\
    , FP_EXT_TOPOLOGICAL_TORSION\
    , FP_MORGAN\
    , FP_TOPOLOGICAL\
    , FP_TOPOLOGICAL_LAYERED_1\
    , FP_TOPOLOGICAL_LAYERED_2\
    , FP_TOPOLOGICAL_TORSION\
    , FP_VECTORFP

FP_ATOM_PAIRS = FP_ATOM_PAIRS
"""
..  todo:: document this
"""

FP_EXT_ATOM_PAIRS = FP_EXT_ATOM_PAIRS
"""
..  todo:: document this
"""

FP_EXT_MORGAN = FP_EXT_MORGAN
"""
..  todo:: document this
"""

FP_EXT_TOPOLOGICAL = FP_EXT_TOPOLOGICAL
"""
..  todo:: document this
"""

FP_EXT_TOPOLOGICAL_LAYERED_1 = FP_EXT_TOPOLOGICAL_LAYERED_1
"""
..  todo:: document this
"""

FP_EXT_TOPOLOGICAL_LAYERED_2 = FP_EXT_TOPOLOGICAL_LAYERED_2
"""
..  todo:: document this
"""

FP_EXT_TOPOLOGICAL_TORSION = FP_EXT_TOPOLOGICAL_TORSION
"""
..  todo:: document this
"""

FP_MORGAN = FP_MORGAN
"""
This is the default selector that is used
when no other option is specified.
"""

FP_TOPOLOGICAL = FP_TOPOLOGICAL
"""
..  todo:: document this
"""

FP_TOPOLOGICAL_LAYERED_1 = FP_TOPOLOGICAL_LAYERED_1
"""
..  todo:: document this
"""

FP_TOPOLOGICAL_LAYERED_2 = FP_TOPOLOGICAL_LAYERED_2
"""
..  todo:: document this
"""

FP_TOPOLOGICAL_TORSION = FP_TOPOLOGICAL_TORSION
"""
..  todo:: document this
"""

FP_VECTORFP = FP_VECTORFP
"""
..  todo:: document this
"""

from molpher.swig_wrappers.core import \
    SC_ALL_BIT\
    , SC_ASYMMETRIC\
    , SC_BRAUN_BLANQUET\
    , SC_COSINE\
    , SC_DICE\
    , SC_KULCZYNSKI\
    , SC_MC_CONNAUGHEY\
    , SC_ON_BIT\
    , SC_RUSSEL\
    , SC_SOKAL\
    , SC_TANIMOTO\
    , SC_TVERSKY_SUBSTRUCTURE\
    , SC_TVERSKY_SUPERSTRUCTURE

SC_ALL_BIT = SC_ALL_BIT
"""
..  todo:: document this
"""

SC_ASYMMETRIC = SC_ASYMMETRIC
"""
..  todo:: document this
"""

SC_BRAUN_BLANQUET = SC_BRAUN_BLANQUET
"""
..  todo:: document this
"""

SC_COSINE = SC_COSINE
"""
..  todo:: document this
"""

SC_DICE = SC_DICE
"""
..  todo:: document this
"""

SC_KULCZYNSKI = SC_KULCZYNSKI
"""
..  todo:: document this
"""

SC_MC_CONNAUGHEY = SC_MC_CONNAUGHEY
"""
..  todo:: document this
"""

SC_ON_BIT = SC_ON_BIT
"""
..  todo:: document this
"""

SC_RUSSEL = SC_RUSSEL
"""
..  todo:: document this
"""

SC_SOKAL = SC_SOKAL
"""
..  todo:: document this
"""

SC_TANIMOTO = SC_TANIMOTO
"""
This is the default selector that is used
when no other option is specified.
"""

SC_TVERSKY_SUBSTRUCTURE = SC_TVERSKY_SUBSTRUCTURE
"""
..  todo:: document this
"""

SC_TVERSKY_SUPERSTRUCTURE = SC_TVERSKY_SUPERSTRUCTURE
"""
..  todo:: document this
"""

from molpher.swig_wrappers.core import \
    OP_ADD_ATOM\
    , OP_ADD_BOND\
    , OP_BOND_CONTRACTION\
    , OP_BOND_REROUTE\
    , OP_INTERLAY_ATOM\
    , OP_MUTATE_ATOM\
    , OP_REMOVE_ATOM\
    , OP_REMOVE_BOND

OP_ADD_ATOM = OP_ADD_ATOM
"""
..  todo:: document this
"""

OP_ADD_BOND = OP_ADD_BOND
"""
..  todo:: document this
"""

OP_BOND_CONTRACTION = OP_BOND_CONTRACTION
"""
..  todo:: document this
"""

OP_BOND_REROUTE = OP_BOND_REROUTE
"""
..  todo:: document this
"""

OP_INTERLAY_ATOM = OP_INTERLAY_ATOM
"""
..  todo:: document this
"""

OP_MUTATE_ATOM = OP_MUTATE_ATOM
"""
..  todo:: document this
"""

OP_REMOVE_ATOM = OP_REMOVE_ATOM
"""
..  todo:: document this
"""

OP_REMOVE_BOND = OP_REMOVE_BOND
"""
..  todo:: document this
"""
