import unittest
from molpher.core import *

class TestPythonAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_source = 'CCO'
        self.test_target = 'O1C=CC=C1'

    def tearDown(self):
        pass

    def testParams(self):
        pass # TODO: test the ExplorationParameters class for completeness

    def testTreeInit(self):
        mol1 = self.test_source
        mol2 = self.test_target
        params = ExplorationParameters(
            source=mol1
            , target=mol2
        )
        params.params = {
            'source' : mol2
            , 'target' : mol1
            , 'operators' : params.params['operators'][:3]
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
        self.assertEqual(tree.params['operators'], params.params['operators'])
        tree.params = {
            'source' : mol1
            , 'target' : mol2
        }
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], mol2)
        _params = ExplorationParameters()
        self.assertEqual(tree.params['operators'], _params.params['operators'])

    def testTreeMethods(self):
        tree = ExplorationTree(params={
            'source' : self.test_source
            , 'target' : self.test_target
        })
        tree.generateMorphs()
        pass
