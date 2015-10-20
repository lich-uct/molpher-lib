
import unittest
from molpher import *

class TestStringMethods(unittest.TestCase):
    
    def setUp(self):
        self.test_files_path = "tests/test_files/"

    def tearDown(self):
        pass


    def testExplorationTreeSnapshotClass(self):
        etreeSnap = ExplorationTreeSnapshot.load(self.test_files_path + "test-template.xml")
        etreeSnap.save(self.test_files_path + "snappitty_snap.snp");

if __name__ == '__main__':
    unittest.main()
