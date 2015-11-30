
import unittest
from pkg_resources import resource_filename
from molpher.swig_wrappers.core import *
import time

class TestMolpherAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_files_path = resource_filename('molpher.swig_wrappers.tests', 'test_files/')

    def tearDown(self):
        pass
        
    def testExplorationParametersClass(self):
        params = ExplorationParameters()
        self.assertFalse(params.valid())
        params.setSourceMol("CCO")
        self.assertTrue(params.valid())
        new_opers = ('ADD_ATOM', 'ADD_BOND', 'BOND_CONTRACTION')
        new_fp = 'EXT_MORGAN'
        new_coef = 'MC_CONNAUGHEY'
        params.setChemOperators(new_opers)
        params.setFingerprint(new_fp)
        params.setSimilarityCoef(new_coef)
        self.assertEqual(params.getChemOperators(), new_opers)
        self.assertEqual(params.getFingerprint(), new_fp)
        self.assertEqual(params.getSimilarityCoef(), new_coef)
        self.assertRaises(RuntimeError, lambda: params.setFingerprint('UNKNOWN'))
        self.assertRaises(RuntimeError, lambda: params.setSimilarityCoef('UNKNOWN'))
        
    def testExplorationTreeSnapshotClass(self):
        etreeSnap = ExplorationTreeSnapshot.load(self.test_files_path + "test-template.xml")
        etreeSnap.save(self.test_files_path + "snappitty_snap.snp")
        
    def testExplorationTreeClass(self):
        params = ExplorationParameters()
        params.setSourceMol("CCO")
        param_tree = ExplorationTree(params)
        snap = param_tree.createSnapshot()
        
        smiles = "OCCO"
        smile_tree = ExplorationTree(smiles)
        snap = smile_tree.createSnapshot()
        snap.save(self.test_files_path + "snappy_snap.snp")
        snap = snap.load(self.test_files_path + "snappy_snap.snp")
        etree = ExplorationTree.createFromSnapshot(snap)
        mol = etree.fetchMol(smiles)
        self.assertEqual(mol.getSMILES(), smiles)
        
    def testExploration(self):
        etreeSnap = ExplorationTreeSnapshot.load(self.test_files_path + "test-template.xml")
        tree = ExplorationTree.createFromSnapshot(etreeSnap)
        tree.setThreadCount(2)

        # set some parameters
        params = tree.getParams()
        new_opers = ('ADD_ATOM', 'ADD_BOND',)
        new_fp = 'EXT_MORGAN'
        new_coef = 'MC_CONNAUGHEY'
        params.setChemOperators(new_opers)
        params.setFingerprint(new_fp)
        params.setSimilarityCoef(new_coef)
        tree.setParams(params)
        
        # find leaves
        print("Searching for leaves...")
        leaves = tree.fetchLeaves()
        self.assertEqual(len(leaves),1)
        self.assertEqual('', leaves[0].getParentSMILES())
        print(leaves)
        print(leaves[0].getSMILES())
        
        # generate morphs
        print("Generating morphs...")
        tree.generateMorphs()
        morphs = tree.getCandidateMorphs()
        self.assertEqual(len(leaves),1)
        for morph in morphs:
            self.assertTrue(morph.getSMILES())
        mask = tree.getCandidateMorphsMask()
        tree.setCandidateMorphsMask(mask)
        self.assertRaises(RuntimeError, lambda: tree.setCandidateMorphsMask([False]))
        self.assertEqual(len(mask), len(morphs))
        print("morphs generated: " + str(len(morphs)))
            
        # sort morphs
        print("Sorting morphs...")
        tree.sortMorphs()
        morphs = tree.getCandidateMorphs()
        previous = None
        for morph in morphs:
            if previous:
                self.assertTrue(morph.getDistToTarget() >= previous.getDistToTarget())
            previous = morph
            
        # filter morphs
        print("Filtering morphs...")
        tree.filterMorphs(FilterMoprhsOper.COUNT | FilterMoprhsOper.WEIGHT | FilterMoprhsOper.PROBABILITY)
        mask = tree.getCandidateMorphsMask()
        print("count - weight - probability filter survivors: " + str(sum(mask)))
        print(mask)
        tree.filterMorphs(FilterMoprhsOper.ALL)
        mask = tree.getCandidateMorphsMask()
        print("all filter survivors: " + str(sum(mask)))
        print(mask)
        self.assertEqual(len(mask), len(morphs))
        
        # extend the tree with the accepted morphs (true in the survivors mask)
        print("Extending tree...")
        tree.extend()
        new_leaves = tree.fetchLeaves()
        self.assertEqual(sum(mask), len(new_leaves))
        for leaf in new_leaves:
            parent_smi = leaf.getParentSMILES()
            self.assertEqual(parent_smi, leaves[0].getSMILES())
            
        self.assertEqual(len(tree.getCandidateMorphsMask()), 0)
        self.assertEqual(len(tree.getCandidateMorphs()), 0)
        
        # remove subtree from tree
        print("Removing subtree...")
        tree.deleteSubtree(new_leaves[0].getSMILES())
        self.assertEqual(len(new_leaves) - 1, len(tree.fetchLeaves()))
        
        # pruning
        print("Pruning tree...")
        tree.prune()
        
    def testTreeOperCallback(self):
        class MyTreeOper(TreeOperation):
            
            def __init__(self):
                super(MyTreeOper, self).__init__()
                self.new_morphs = None
            
            def __call__(self):
                tree = self.getTree()
                tree.generateMorphs()
                self.new_morphs = [ x.getSMILES() for x in tree.getCandidateMorphs()]
                
        tree = ExplorationTree("CCO")
        oper = MyTreeOper()
        self.assertEqual(oper.new_morphs, None)
        tree.runOperation(oper)
        new_morphs = oper.new_morphs
        tree_morphs = [ x.getSMILES() for x in tree.getCandidateMorphs()]
        print(new_morphs)
        print(tree_morphs)
        self.assertEqual(tree_morphs, new_morphs)
        
    def testTreeTraversalCallback(self):
        class MyCallback(TraverseCallback):

            def __init__(self, tree):
                super(self.__class__, self).__init__()
                self.tree = tree
                params = self.tree.getParams()
                self.source = params.getSourceMol()
                self.target = params.getTargetMol()
                self.already_printed_stuff = False
            
            def processMorph(self, morph):
                time.sleep(0.020)
                morph.getSMILES()
                if morph.getHistoricDescendants() and not self.already_printed_stuff:
                    print("Found root.")
                    print("Historic descendants:")
                    print(morph.getHistoricDescendants())
                    descendants = morph.getDescendants()
                    print("Descendants:")
                    print(descendants)
                    for desc_smile in descendants:
                        desc_mol = self.tree.fetchMol(desc_smile)
                        print("Descendant {0} from source {1} to target {2}:".format(desc_smile, self.source.getSMILES(), self.target.getSMILES()))
                        print("Weight: {0}".format(desc_mol.getMolecularWeight()))
                        print("SAScore: {0}".format(desc_mol.getSAScore()))
                        print("Distance to target: {0}".format(desc_mol.getDistToTarget()))

                    self.already_printed_stuff = True
                
                
        tree = ExplorationTree("CCO")
        tree.setThreadCount(2)
        tree.generateMorphs()
        tree.filterMorphs(FilterMoprhsOper.WEIGHT)
        tree.extend()
        callback = MyCallback(tree)
        traverse = TraverseOper(tree, callback)
        count = 0
        while count < 30:
            traverse()
            count += 1

if __name__ == '__main__':
    unittest.main()
