from molpher.algorithms.functions import find_path
from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

class NitorgenFilter(TreeOperation):

    def __call__(self):
        new_mask = [ 'N' in x.smiles for x in self.tree.candidates ]
        self.tree.candidates_mask = new_mask


iteration = [
    GenerateMorphsOper()
    , SortMorphsOper()
    , FilterMorphsOper() # the default filter
    , CleanMorphsOper() # discards morphs that were previously filtered out
    , NitorgenFilter() # customized filter
    , ExtendTreeOper()
    , PruneTreeOper()
]

tree = ETree.create(source=cocaine, target=procaine)
counter = 0
while not tree.path_found:
    counter+=1
    print("Iteration", counter)
    for oper in iteration:
        tree.runOperation(oper)

print("Path found: ", find_path(tree))


