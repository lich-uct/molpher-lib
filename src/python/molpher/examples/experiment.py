"""
Complete morphing example create from the tutorial code.

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

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

class MyFilterMorphs(TreeOperation):

    def __init__(self):
        super(MyFilterMorphs, self).__init__()

    def __call__(self):
        mask = [False for x in self.tree.candidates_mask]
        mask[0] = True
        mask[1] = True
        mask[2] = True
        self.tree.candidates_mask = mask

def main():
    iteration = [
        GenerateMorphsOper()
        , SortMorphsOper()
        , MyFilterMorphs()
        , ExtendTreeOper()
        , PruneTreeOper()
    ]

    tree = ETree.create(source=cocaine, target=procaine)
    counter = 0
    while not tree.path_found:
        for oper in iteration:
            tree.runOperation(oper)
        counter+=1
        print("Iteration", counter)
        print(
            sorted(
            [
                (x.getSMILES(), x.getDistToTarget())
                for x in tree.leaves
            ], key=lambda x : x[1]
            )
        )

    path = []
    current = tree.fetchMol(tree.params['target'])
    path.append(current)
    while current != '':
        current = current.getParentSMILES()
        if current:
            current = tree.fetchMol(current)
            path.append(current)

    path.reverse()
    print("Path found: ")
    for mol in path:
        print(mol.getSMILES(), mol.getDistToTarget())


if __name__ == "__main__":
    exit(main())