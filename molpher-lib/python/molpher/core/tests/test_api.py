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

    def testMorphing(self):
        def callback(morph):
            if morph.getItersWithoutDistImprovement() > 3:
                print(morph.getSMILES(), morph.getItersWithoutDistImprovement())

        class RunIteration(TreeOperation):

            def __init__(self, tree):
                super(RunIteration, self).__init__()
                self.tree = tree

            def __call__(self):
                print('Iteration: ', self.tree.getGenerationCount() + 1)
                self.tree.generateMorphs()
                self.tree.extend()
                self.tree.traverse(callback)

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


