"""
Module with example code from the tutorial.

"""

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

# TODO: make some of the stuff from this script part of the test suite

def main():
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    tree = ETree.create(source=cocaine, target=procaine) # initialize a tree that searches for a path from cocaine to procaine

    # print the smiles of the source and target molecule
    print('Source: ', tree.params['source'])
    print('Target: ', tree.params['target'])

    # change selected parameters using a dictionary
    print(tree.params)
    tree.params = {
        'non_producing_survive' : 2,
        'weight_max' : 500.0
    }
    print(tree.params)

    print('\n#Generating and Manipulating Morphs')

    print(tree.leaves) # show the current leaves of the tree (only the source so far)
    print(tree.leaves[0].smiles)
    tree.generateMorphs() # generate new morphs
    print(tree.candidates)
    print(len(tree.candidates))

    print()

    # get the first morph in the candidate list
    candidate = tree.candidates[0]
    # print distance to target
    print(tree.candidates[0].dist_to_target)
    # set new distance to target
    candidate.dist_to_target = 0.5
    # look in the list of candidates and print new distance
    print(tree.candidates[0].dist_to_target)

    print()

    # make a copy of our molecule
    candidate_copy = candidate.copy()
    # set a new distance for the copy and verify that the original was not affected
    print(candidate_copy.dist_to_target)
    candidate_copy.dist_to_target = 0.7
    print(candidate_copy.dist_to_target)
    print(candidate.dist_to_target)
    print(tree.candidates[0].dist_to_target)

    print('\n#Sorting and Filtering Morphs')

    # sort the candidates in the tree according to their distance from target
    tree.sortMorphs()

    # show results
    print(tree.candidates_mask)
    print(
        [
            (x.smiles, x.dist_to_target)
            for x in tree.candidates
        ]
    )

    print()

    # print the current candidates mask (all positions are on by default)
    print(tree.candidates_mask)

    # accept only the first three morphs in the sorted list (those with the lowest distance to target)
    mask = [False for x in tree.candidates_mask]
    mask[0] = True
    mask[1] = True
    mask[2] = True

    # save the new mask to the tree
    tree.candidates_mask = mask

    # show results
    print(tree.candidates_mask)
    print(
        [
            (x.smiles, x.dist_to_target)
            for idx,x in enumerate(tree.candidates)
            if tree.candidates_mask[idx] # get accepted molecules only
        ]
    )

    print('\n#Extending and Pruning')

    # get the number of generations before
    print(tree.generation_count)

    tree.extend() # connect the accepted morphs to the tree as new leaves
    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda item : item[1]
        )
    )

    # get the number of generations after
    print(tree.generation_count)

    # check if a path was found
    print(tree.path_found)

    # run the pruning operation on the updated tree
    tree.prune()

    print('\n#Operations')

    class MyFilterMorphs(TreeOperation):
        """
        A custom tree operation that accepts
        only the first three morphs
        (those with the lowest distance to target).

        """

        def __call__(self):
            """
            This method is called automatically by the tree.
            The tree this operation is being run on is accessible
            from `self.tree`.

            """

            mask = [False for x in self.tree.candidates_mask]
            mask[0] = True
            mask[1] = True
            mask[2] = True
            self.tree.candidates_mask = mask

    tree = ETree.create(source=cocaine, target=procaine) # create the tree

    # this list of tree operations defines one iteration
    iteration = [
        GenerateMorphsOper()
        , SortMorphsOper()
        , MyFilterMorphs()
        , ExtendTreeOper()
        , PruneTreeOper()
    ]

    # apply the operations in the list one by one
    for oper in iteration:
        tree.runOperation(oper)

    # observe the results
    print(tree.generation_count)
    print(tree.path_found)
    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )

    print('\n#Traversing the Tree')

    class MyCallback(TraverseCallback):
        """
        This callback just prints some information
        about the molecules in the tree.

        """

        def __call__(self, morph):
            """
            Method called on each morph in the tree
            -- starting from root to leaves.

            """

            if not morph.getParentSMILES():
                print("# Root #")
            else:
                print('# Morph #')
                print('Parent:', morph.getParentSMILES())
            print('SMILES: ', morph.getSMILES())
            print('Descendents: ', morph.getDescendants())

    callback = MyCallback() # initialize a callback
    traverse = TraverseOper(callback) # attach it to a tree traversal operation
    tree.runOperation(traverse) # run the operation

    print()

    def process(morph):
        """
        Prints some information
        about the molecules in the tree.

        """

        if not morph.getParentSMILES():
            print("# Root #")
        else:
            print('# Morph #')
            print('Parent:', morph.getParentSMILES())
        print('SMILES: ', morph.getSMILES())
        print('Descendents: ', morph.getDescendants())

    tree.traverse(process) # use the traverse method to run the callback function

    print('\n#Tree Snapshots')

    template_file = 'cocaine-procaine-template.xml'
    import os
    template_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), template_file)

    # create a tree from the template file
    tree = ETree.create(template_file)
    print(tree.params)

    # apply the tree operations
    for oper in iteration:
        tree.runOperation(oper)

    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )

    # save the tree in a snapshot file
    tree.save('snapshot.xml')

    new_tree = ETree.create('snapshot.xml') # create a new tree from the saved snapshot
    print(new_tree.params)
    print(
        sorted( # grab the leaves in the created tree (these should be the same as those in the original tree)
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in new_tree.leaves
        ], key=lambda x : x[1]
        )
    )

if __name__ == "__main__":
    exit(main())