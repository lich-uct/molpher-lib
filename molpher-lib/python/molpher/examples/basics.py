from molpher.core import ExplorationTree as ETree

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