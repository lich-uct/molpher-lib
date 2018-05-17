"""
Complete morphing example created from the tutorial code.

"""

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

class MyFilterMorphs(TreeOperation):

    def __call__(self):

        self.tree.candidates_mask = [
            True if idx < 10 and self.tree.candidates[idx].sascore < 6
            else False
            for idx, x in enumerate(self.tree.candidates_mask)
        ]

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

    path = tree.fetchPathTo(tree.params['target'])
    print("Path found: ")
    for mol in path:
        print(mol.getSMILES(), mol.getDistToTarget())


if __name__ == "__main__":
    exit(main())