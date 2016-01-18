"""
Module with example code from the tutorial.

"""

from molpher.core import ExplorationTree as ETree
from molpher.core.operations import *

# TODO: make some of the stuff from this script part of the test suite

def main():
    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    tree = ETree(source=cocaine, target=procaine) # initialize a tree that searches for a path from cocaine to procaine

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

    print('#Generating and Manipulating Morphs')

    print(tree.leaves) # show the current leaves of the tree (only the source so far)
    print(tree.leaves[0].getSMILES())
    tree.generateMorphs() # generate new morphs
    print(tree.candidates)
    print(len(tree.candidates))

    print()

    candidate = tree.candidates[0] # get the first morph in the candidate list
    print(candidate.isBound()) # print if it is bound to a tree
    print(tree.candidates[0].getDistToTarget()) # print distance to target
    candidate.setDistToTarget(0.5) # set new distance to target
    print(tree.candidates[0].getDistToTarget()) # look in the list of candidates and print new distance

    print()

    # get a copy of the candidate molecule and verify that the changes do not propagate into the tree
    candidate_copy = candidate.copy()
    print(candidate_copy.isBound())
    print(candidate_copy.getDistToTarget())
    candidate_copy.setDistToTarget(0.7)
    print(candidate_copy.getDistToTarget())
    print(candidate.getDistToTarget())
    print(tree.candidates[0].getDistToTarget())

    print('#Sorting and Filtering Morphs')

    tree.sortMorphs() # sort the candidates in the tree according to their distance from target
    print(tree.candidates_mask) # print the current mask
    print(len(tree.candidates_mask))

    # accept only the first three morphs (those with the lowest distance to target)
    mask = [False for x in tree.candidates_mask]
    mask[0] = True
    mask[1] = True
    mask[2] = True
    tree.candidates_mask = mask # save the new mask to the tree

    # verify
    print(tree.candidates_mask)
    print(
        [
            (x.getSMILES(), x.getDistToTarget())
            for idx,x in enumerate(tree.candidates)
            if tree.candidates_mask[idx] # only fetch accepted molecules
        ]
    )

    print('#Extending and Pruning')

    tree.extend() # connect the accepted morphs to the tree as new leaves
    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )
    print(tree.generation_count) # get the number of generations
    print(tree.path_found) # check if a path was found
    tree.prune() # run the pruning operation on the updated tree

    print('#Operations')

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

    tree = ETree(source=cocaine, target=procaine) # create the tree

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

    print('Traversing the Tree')

    class MyCallback(TraverseCallback):
        """
        This callback just prints some information
        about the molecules in the tree.

        """

        def processMorph(self, morph):
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

    print('Tree Snapshots')

    template_file = 'cocaine-procaine-template.xml'
    import os
    template_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), template_file)

    tree = ETree.createFromSnapshot(template_file) # create a tree from the template file
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

    tree.saveSnapshot('snapshot.xml') # save the tree in a snapshot file

    new_tree = ETree.createFromSnapshot('snapshot.xml') # create a new tree from the saved snapshot
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