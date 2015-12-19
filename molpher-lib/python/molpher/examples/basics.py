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
print(tree.candidates_mask)
print(len(tree.candidates))
print(tree.leaves[0].isBound())
print(len(tree.candidates_mask))