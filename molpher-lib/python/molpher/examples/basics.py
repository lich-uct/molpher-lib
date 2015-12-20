from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'
cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'

tree = ETree(source=cocaine, target=procaine)

print('Source: ', tree.params['source'])
print('Target: ', tree.params['target'])

print(tree.params)
tree.params = {
    'non_producing_survive' : 2,
    'weight_max' : 500.0
}
print(tree.params)

print()

print(tree.leaves)
print(tree.leaves[0].getSMILES())
tree.generateMorphs()
print(tree.candidates)
print(len(tree.candidates))

print()

candidate = tree.candidates[0]
print(candidate.isBound())
print(tree.candidates[0].getDistToTarget())
candidate.setDistToTarget(0.5)
print(tree.candidates[0].getDistToTarget())

print()

candidate_copy = candidate.copy()
print(candidate_copy.isBound())
print(candidate_copy.getDistToTarget())
candidate_copy.setDistToTarget(0.7)
print(candidate_copy.getDistToTarget())
print(candidate.getDistToTarget())

print()

print(tree.candidates_mask)
print(len(tree.candidates_mask))
mask = [False for x in tree.candidates_mask]
mask[0] = True
mask[1] = True
mask[2] = True
tree.candidates_mask = mask
print(tree.candidates_mask)

print()

print(sorted([x.getSMILES() for idx,x in enumerate(tree.candidates) if tree.candidates_mask[idx]]))
tree.extend()
print(sorted([x.getSMILES() for x in tree.leaves]))
print(tree.generation_count)
print(tree.path_found)
tree.prune()

print()

class MyFilterMorphs(TreeOperation):

    def __init__(self):
        super(MyFilterMorphs, self).__init__()

    def __call__(self):
        tree = self.getTree()
        mask = [False for x in tree.candidates_mask]
        mask[0] = True
        mask[1] = True
        mask[2] = True
        tree.candidates_mask = mask

    def getTree(self):
        tree = super(MyFilterMorphs, self).getTree()
        if tree:
            tree.__class__ = ETree # 'cast' the wrapped class to the pretty Python proxy class
        return tree

tree = ETree(source=cocaine, target=procaine)

iteration = [
    GenerateMorphsOper()
    , SortMorphsOper()
    , MyFilterMorphs()
    , ExtendTreeOper()
    , PruneTreeOper()
]

for oper in iteration:
    tree.runOperation(oper)

print(tree.generation_count)
print([(x.getSMILES(), x.getDistToTarget()) for x in tree.leaves])

# TODO: make some of the stuff from this script part of the test suite