import unittest
from molpher.core import *
from molpher.swig_wrappers.core import TreeOperation

class TestPythonAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_source = 'CCO'
        self.test_target = 'O1C=CC=C1'

    def tearDown(self):
        pass

    def testParams(self):
        params = ExplorationParameters(
            source=self.test_source
            , target=self.test_target
        )
        params2 = ExplorationParameters(parameters=params)
        self.assertEqual(params.source.getSMILES(), params2.source.getSMILES())
        self.assertEqual(params.target.getSMILES(), params2.target.getSMILES())
        self.assertEqual(params.is_valid == True, params2.is_valid == True)

    def testTreeInit(self):
        mol1 = self.test_source
        mol2 = self.test_target
        params = ExplorationParameters(
            source=mol1
            , target=mol2
        )
        params.param_dict = {
            'source' : mol2
            , 'target' : mol1
            , 'operators' : params.param_dict['operators'][:3]
        }
        self.assertRaises(RuntimeError, lambda : ExplorationTree())
        tree = ExplorationTree(source=mol1, target=mol2)
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], mol2)
        def func():
            tree.params = {'source' : ''}
        self.assertRaises(RuntimeError, func)
        tree.params = params
        self.assertEqual(tree.params['source'], mol2)
        self.assertEqual(tree.params['target'], mol1)
        self.assertEqual(tree.params['operators'], params.param_dict['operators'])
        tree.params = {
            'source' : mol1
            , 'target' : mol2
        }
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], mol2)
        self.assertEqual(tree.params['operators'], params.param_dict['operators'])

        leaf = tree.leaves[0]
        self.assertTrue(leaf.isBound())
        leaf.setDistToTarget(0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)

        leaf_copy = leaf.copy()
        self.assertFalse(leaf_copy.isBound())
        self.assertEqual(leaf_copy.getDistToTarget(), 0.5)
        leaf_copy.setDistToTarget(0.7)
        self.assertEqual(leaf.getDistToTarget(), 0.5)
        self.assertEqual(leaf_copy.getDistToTarget(), 0.7)

    def testMorphing(self):
        def callback(morph):
            callback.morphs_in_tree += 1
            self.assertTrue(morph.isBound())
            if morph.getItersWithoutDistImprovement() > 3:
                print(morph.getSMILES(), morph.getItersWithoutDistImprovement())
            if not callback.closest_mol:
                callback.closest_mol = morph.copy()
            if callback.closest_mol.getDistToTarget() > morph.getDistToTarget():
                callback.closest_mol = morph.copy()
        callback.morphs_in_tree = 0
        callback.closest_mol = None

        class RunIteration(TreeOperation):

            def __init__(self, tree):
                super(RunIteration, self).__init__()
                self._tree = tree

            def __call__(self):
                print('Iteration: ', self._tree.getGenerationCount() + 1)
                self._tree.generateMorphs()
                self._tree.filterMorphs()
                self._tree.extend()
                self._tree.prune()
                callback.morphs_in_tree = 0
                self._tree.traverse(callback)
                print('Number of morphs in the tree: ', callback.morphs_in_tree)
                print('Closest molecule to target: {0} -- {1}'.format(
                    callback.closest_mol.getSMILES()
                    , callback.closest_mol.getDistToTarget()
                ))

            def getTree(self):
                return self._tree

            def setTree(self, tree):
                self._tree = tree

        tree = ExplorationTree(params={
            'source' : self.test_source
            , 'target' : self.test_target
        })

        iterate = RunIteration(tree)
        counter = 0
        while counter < 5:
            iterate()
            counter += 1

    def testOperations(self):
        # TODO: implement
        pass

    def testSnapshots(self):
        # TODO: implement
        pass


