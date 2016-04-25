import os
import unittest

from molpher.core import ExplorationTree
from molpher.core import MolpherMol
from molpher.core.operations import TreeOperation
from molpher.core.selectors import *
from molpher.core import ExplorationData

class TestPythonAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_source = 'CCO'
        self.test_target = 'C1=COC=C1'
        self.test_dir = os.path.abspath('test_files/')
        self.test_template_path = os.path.join(self.test_dir, 'test-template.xml')

    def tearDown(self):
        pass

    def testMolpherMol(self):
        mol = MolpherMol(self.test_target)
        mol.smiles = 'CCC'
        self.assertTrue('CCC')
        mol.historic_descendents = ('CCC', 'CCCC')
        self.assertEqual(('CCC', 'CCCC'), mol.historic_descendents)

        copy = mol.copy()
        copy.sascore = 0.54
        self.assertEqual(0.54, copy.sascore)
        self.assertEqual(0, mol.sascore)

    def testExplorationData(self):
        params = ExplorationData(
            source=self.test_source
            , target=self.test_target
        )

        params.operators = set((OP_ADD_BOND, 'OP_REMOVE_BOND'))
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
        # self.assertEqual(tree, leaf.tree) # FIXME: add a reliable operator for comparison between trees
        leaf.setDistToTarget(0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)

        leaf_copy = tree.leaves[0].copy()
        self.assertFalse(tree.hasMol(leaf_copy))
        self.assertEqual(leaf_copy.getDistToTarget(), 0.5)
        leaf_copy.setDistToTarget(0.7)
        self.assertEqual(leaf.getDistToTarget(), 0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)
        self.assertEqual(leaf_copy.getDistToTarget(), 0.7)

    def testOperations(self):
        # TODO: implement
        pass

    def testMorphing(self):
        def callback(morph):
            callback.morphs_in_tree += 1
            self.assertTrue(morph)
            self.assertTrue(morph.tree)
            if morph.getItersWithoutDistImprovement() > 3:
                print('Callback output:')
                print(morph.getSMILES(), morph.getItersWithoutDistImprovement(), morph.getDistToTarget())
            if not callback.closest_mol:
                callback.closest_mol = morph.copy()
            current_dist = morph.getDistToTarget()
            min_dist = callback.closest_mol.getDistToTarget()
            if min_dist > current_dist:
                callback.closest_mol = morph.copy()

            # FIXME: code below doesn't work
                # (causes a SEGFAULT when the memory is accessed,
                # probably because the shared pointer is deleted
                # from the stack in the SWIG generated code when the callback is done)

            # if not callback.closest_mol:
            #     callback.closest_mol = morph
            # current_dist = morph.getDistToTarget()
            # min_dist = callback.closest_mol.getDistToTarget()
            # if min_dist > current_dist:
            #     callback.closest_mol = morph
        callback.morphs_in_tree = 0
        callback.closest_mol = None

        class MorphingIteration(TreeOperation):

            parent = self

            def __init__(self, tree):
                super(MorphingIteration, self).__init__()
                self._tree = tree

            def __call__(self):
                print('Iteration: ', self._tree.getGenerationCount() + 1)
                self._tree.generateMorphs()
                for mol in self._tree.candidates:
                    self.parent.assertEqual(None, mol.tree)
                self._tree.filterMorphs()
                self._tree.extend()
                self._tree.prune()
                callback.morphs_in_tree = 0
                self._tree.traverse(callback)
                print('Number of morphs in the tree: ', callback.morphs_in_tree)
                print('Closest molecule to target: {0} -- distance: {1}'.format(
                    callback.closest_mol.getSMILES()
                    , callback.closest_mol.getDistToTarget()
                ))

            def getTree(self):
                return self._tree

            def setTree(self, tree):
                self._tree = tree

        tree = ExplorationTree.create(params={
            'source' : self.test_source
            , 'target' : self.test_target
        })

        iterate = MorphingIteration(tree)
        counter = 0
        while counter < 5:
            iterate()
            counter += 1

        child = tree.leaves[0]
        self.assertTrue(child.tree)
        self.assertTrue(tree.hasMol(child))
        parent = child.getParentSMILES()
        tree.deleteSubtree(parent)
        self.assertFalse(tree.hasMol(parent))
        self.assertFalse(tree.hasMol(child))
        self.assertEqual(None, child.tree)
        self.assertEqual(parent, child.getParentSMILES())