import os
import unittest

from molpher.core.ExplorationTree import ExplorationTree
from molpher.core.selectors import *
from molpher.core.ExplorationData import ExplorationData

class TestPythonAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_source = 'CCO'
        self.test_target = 'C1=COC=C1'
        self.test_dir = os.path.abspath('test_files/')
        self.test_template_path = os.path.join(self.test_dir, 'test-template.xml')

    def tearDown(self):
        pass

    def test_MolpherMol(self):
        # TODO: implement
        pass

    def testParams(self):
        params = ExplorationData(
            source=self.test_source
            , target=self.test_target
        )

        params.operators = {OP_ADD_BOND, 'OP_REMOVE_BOND'}
        self.assertEqual(params.operators, ('OP_ADD_BOND', 'OP_REMOVE_BOND'))

        params.fingerprint = FP_EXT_ATOM_PAIRS
        self.assertEqual(params.fingerprint, 'FP_EXT_ATOM_PAIRS')
        params.fingerprint = 'FP_TOPOLOGICAL_LAYERED_2'
        self.assertEqual(params.fingerprint, 'FP_TOPOLOGICAL_LAYERED_2')
        params.similarity = 'SC_COSINE'
        self.assertEqual(params.similarity, 'SC_COSINE')
        params.similarity = SC_KULCZYNSKI
        self.assertEqual(params.similarity, 'SC_KULCZYNSKI')

        self.assertEqual(params.source.smiles, self.test_source)
        self.assertEqual(params.target.smiles, self.test_target)

        params.param_dict = {
            'target' : self.test_source
            , 'operators' : params.param_dict['operators'][:1]
        }
        self.assertEqual(params.target.smiles, self.test_source)
        self.assertEqual(params.operators, ('OP_ADD_BOND',))

        params_from_temp = ExplorationData.load(self.test_template_path)
        self.assertRaises(RuntimeError, lambda : ExplorationData.load('not/a/valid/path'))
        self.assertEqual(777, params_from_temp.far_produce)

    def testTree(self):
        mol1 = self.test_source
        mol2 = self.test_target
        params_dict = {
            'source' : mol1
            , 'target' : mol2
            , 'operators' : (OP_ADD_BOND, OP_REMOVE_BOND, OP_MUTATE_ATOM)
        }
        params = ExplorationData(**params_dict)


        self.assertRaises(AttributeError, lambda : ExplorationTree())
        tree_from_dict = ExplorationTree.create(params=params_dict)
        tree_from_params = ExplorationTree.create(params=params)
        tree_from_SMILES = ExplorationTree.create(source=mol1, target=mol2)
        def test_tree(tree):
            self.assertEqual(tree.params['source'], mol1)
            self.assertEqual(tree.params['target'], mol2)

        test_tree(tree_from_dict)
        test_tree(tree_from_params)
        test_tree(tree_from_SMILES)

        tree = tree_from_params

        def func():
            tree.params = {'source' : ''}
        self.assertRaises(RuntimeError, func)

        tree.params = {
            'source' : mol2 # if we try to set source for non-empty tree, it shouldn't be changed
            , 'target' : 'C'
        }
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], 'C')
        self.assertEqual(tree.params['operators'], params.param_dict['operators']) # we should still have the same opers set

        tree.params = params # assign the original parameters back
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], mol2)
        self.assertEqual(tree.params['operators'], params.param_dict['operators'])

        leaf = tree.leaves[0]
        self.assertTrue(tree.hasMol(leaf))
        # self.assertEqual(tree, leaf.tree) # FIXME: add a reliable operator for tree comparisons
        leaf.setDistToTarget(0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)

        leaf_copy = tree.leaves[0].copy()
        self.assertFalse(tree.hasMol(leaf_copy))
        self.assertEqual(leaf_copy.getDistToTarget(), 0.5)
        leaf_copy.setDistToTarget(0.7)
        self.assertEqual(leaf.getDistToTarget(), 0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)
        self.assertEqual(leaf_copy.getDistToTarget(), 0.7)

    # def testMorphing(self):
    #     def callback(morph):
    #         callback.morphs_in_tree += 1
    #         self.assertTrue(morph.isBound())
    #         if morph.getItersWithoutDistImprovement() > 3:
    #             print(morph.getSMILES(), morph.getItersWithoutDistImprovement())
    #         if not callback.closest_mol:
    #             callback.closest_mol = morph.copy()
    #         if callback.closest_mol.getDistToTarget() > morph.getDistToTarget():
    #             callback.closest_mol = morph.copy()
    #     callback.morphs_in_tree = 0
    #     callback.closest_mol = None
    #
    #     class RunIteration(TreeOperation):
    #
    #         def __init__(self, tree):
    #             super(RunIteration, self).__init__()
    #             self._tree = tree
    #
    #         def __call__(self):
    #             print('Iteration: ', self._tree.getGenerationCount() + 1)
    #             self._tree.generateMorphs()
    #             self._tree.filterMorphs()
    #             self._tree.extend()
    #             self._tree.prune()
    #             callback.morphs_in_tree = 0
    #             self._tree.traverse(callback)
    #             print('Number of morphs in the tree: ', callback.morphs_in_tree)
    #             print('Closest molecule to target: {0} -- {1}'.format(
    #                 callback.closest_mol.getSMILES()
    #                 , callback.closest_mol.getDistToTarget()
    #             ))
    #
    #         def getTree(self):
    #             return self._tree
    #
    #         def setTree(self, tree):
    #             self._tree = tree
    #
    #     tree = ExplorationTree(params={
    #         'source' : self.test_source
    #         , 'target' : self.test_target
    #     })
    #
    #     iterate = RunIteration(tree)
    #     counter = 0
    #     while counter < 5:
    #         iterate()
    #         counter += 1
    #
    # def testOperations(self):
    #     # TODO: implement
    #     pass
    #
    # def testSnapshots(self):
    #     # TODO: implement
    #     pass


