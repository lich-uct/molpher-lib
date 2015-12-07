import unittest
from molpher.core import *

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
        self.assertRaises(NotImplementedError, lambda : ExplorationTree(source=mol1, target=mol2))
        self.assertRaises(RuntimeError, lambda : ExplorationTree())
        tree = ExplorationTree(source=mol1)
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], 'C')
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

    def testTreeMethods(self):
        tree = ExplorationTree(params={
            'source' : self.test_source
            , 'target' : self.test_target
        })
        tree.generateMorphs()
        pass
