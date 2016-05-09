"""
Implementation of the bidirectional search algorithm.

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

from molpher.core.ExplorationTree import ExplorationTree as ETree
from molpher.core.operations import *

def timeit(func):
    milliseconds = 1000 * time.clock()
    func()
    return 1000 * time.clock() - milliseconds

class FindClosest:

    def __init__(self):
        self.closest = None

    def __call__(self, morph):
        if not self.closest:
            self.closest = morph.copy()
            return
        current_dist = self.closest.getDistToTarget()
        morph_dist = morph.getDistToTarget()
        if morph_dist < current_dist:
            self.closest = morph.copy()

class BidirectionalPathFinder:

    def __init__(self, source, target):
        options = {
            'fingerprint' : 'ATOM_PAIRS'
        }
        self.source_target = ETree.create(source=source, target=target)
        self.target_source = ETree.create(source=target, target=source)
        self.source_target.params = options
        self.target_source.params = options
        self.source_target_min = FindClosest()
        self.target_source_min = FindClosest()
        self.ITERATION = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , ExtendTreeOper()
            , PruneTreeOper()
        ]
        self.path = []

    def _find_path(self, tree, connecting_mol):
        path = []
        current = tree.fetchMol(connecting_mol)
        path.append(current.getSMILES())
        while current != '':
            current = current.getParentSMILES()
            if current:
                current = tree.fetchMol(current)
                path.append(current.getSMILES())
        path.reverse()
        return path

    def __call__(self):
        counter = 0
        connecting_molecule = None
        while True:
            counter+=1
            print('Iteration {0}:'.format(counter))
            for oper in self.ITERATION:
                print('Execution times ({0}):'.format(type(oper).__name__))

                source_target_time = timeit(lambda : self.source_target.runOperation(oper))
                print('\tsource -> target: {0}'.format(source_target_time))
                target_source_time = timeit(lambda : self.target_source.runOperation(oper))
                print('\ttarget -> source: {0}'.format(target_source_time))

                print('\ttotal time: {0}'.format(source_target_time + target_source_time))

            print('Traversal times:')

            source_target_time = timeit(lambda : self.source_target.traverse(self.source_target_min))
            print('\tsource -> target: {0}'.format(source_target_time))
            target_source_time = timeit(lambda : self.target_source.traverse(self.target_source_min))
            print('\ttarget -> source: {0}'.format(target_source_time))

            print('\ttotal time: {0}'.format(source_target_time + target_source_time))

            print('Current Targets:')
            print('\tsource to target:', self.source_target.params['target'])
            print('\ttarget to source:', self.target_source.params['target'])

            print('Current Minima:')
            print('\tsource to target:', self.source_target_min.closest.getSMILES(), self.source_target_min.closest.getDistToTarget())
            print('\ttarget to source:', self.target_source_min.closest.getSMILES(), self.target_source_min.closest.getDistToTarget())

            self.source_target.params = {
                'target' : self.target_source_min.closest.getSMILES()
            }
            self.target_source.params = {
                'target' : self.source_target_min.closest.getSMILES()
            }

            print('New Targets:')
            print('\tsource to target:', self.source_target.params['target'])
            print('\ttarget to source:', self.target_source.params['target'])

            if self.source_target.path_found:
                print('Path Found in tree going from source to target:')
                connecting_molecule = self.source_target.params['target']
                print('Connecting molecule:', connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                break
            if self.target_source.path_found:
                print('Path Found in tree going from target to source:')
                connecting_molecule = self.target_source.params['target']
                print('Connecting molecule:', connecting_molecule)
                assert self.target_source.hasMol(connecting_molecule)
                assert self.source_target.hasMol(connecting_molecule)
                break

        source_target_path = self._find_path(self.source_target, connecting_molecule)
        target_source_path = self._find_path(self.target_source, connecting_molecule)
        assert source_target_path.pop(-1) == connecting_molecule
        target_source_path.reverse()
        source_target_path.extend(target_source_path)
        self.path = source_target_path

def main():
    milliseconds_now = 1000 * time.clock()
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    pathfinder = BidirectionalPathFinder(cocaine, procaine)
    pathfinder()
    print(pathfinder.path)
    print('Total Execution Time: {0}'.format(1000 * time.clock() - milliseconds_now))

if __name__ == "__main__":
    exit(main())