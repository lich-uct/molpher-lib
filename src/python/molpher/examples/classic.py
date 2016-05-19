"""
Implementation of the classical search algorithm.

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

import time

import gc

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds

class ClassicPathFinder:

    def __init__(self, source, target):
        options = {
            'fingerprint' : 'ATOM_PAIRS'
        }
        self.tree = ETree.create(source=source, target=target)
        self.ITERATION = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

    @property
    def path(self):
        current = self.tree.fetchMol(self.tree.params['target'])
        path = list()
        path.append(current.smiles)
        while current != '':
            current = current.parent_smiles
            if current:
                current = self.tree.fetchMol(current)
                path.append(current.smiles)

        path.reverse()
        return path

    def __call__(self):
        counter = 0
        connecting_molecule = None
        while not self.tree.path_found:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                print('Execution time ({0}):'.format(type(oper).__name__))

                time_elapsed = timeit(lambda : self.tree.runOperation(oper))
                print('\telapsed time: {0}'.format(time_elapsed))

def main():
    milliseconds_now = 1000 * time.clock()
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    pathfinder = ClassicPathFinder(cocaine, procaine)
    pathfinder()
    print(pathfinder.path)
    print('Total Execution Time: {0}'.format(1000 * time.clock() - milliseconds_now))

if __name__ == "__main__":
    exit(main())