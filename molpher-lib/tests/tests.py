
import unittest
from molpher import *

class TestMolpherAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_files_path = "tests/test_files/"

    def tearDown(self):
        pass

    def testMolpherMolClass(self):
        mol = MolpherMol("CCO")
        self.assertEqual(mol.getSMILES(), "CCO")
        
    def testExplorationParametersClass(self):
        params = ExplorationParameters()
        self.assertFalse(params.valid())
        params.setSourceMol("CCO")
        self.assertTrue(params.valid())
        
    def testExplorationTreeSnapshotClass(self):
        etreeSnap = ExplorationTreeSnapshot.load(self.test_files_path + "test-template.xml")
        etreeSnap.save(self.test_files_path + "snappitty_snap.snp");
        
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
        
        # find leaves
        leaves = tree.fetchLeaves()
        self.assertEqual(len(leaves),1)
        self.assertEqual('', leaves[0].getParentSMILES())
        
        # generate morphs
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
        tree.sortMorphs()
        morphs = tree.getCandidateMorphs()
        previous = None
        for morph in morphs:
            if previous:
                self.assertTrue(morph.getDistToTarget() >= previous.getDistToTarget())
            previous = morph
            
        # filter morphs
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
        tree.extend()
        new_leaves = tree.fetchLeaves()
        self.assertEqual(sum(mask), len(new_leaves))
        for leaf in new_leaves:
            parent_smi = leaf.getParentSMILES()
            self.assertEqual(parent_smi, leaves[0].getSMILES())

if __name__ == '__main__':
    unittest.main()
