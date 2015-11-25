import unittest

import shutil

from molpher.swig_wrappers.core import *
from pkg_resources import resource_filename

class TestMorphing(unittest.TestCase):

    def setUp(self):
        self.test_files_path = resource_filename('molpher.swig_wrappers', 'test_files/')

    def tearDown(self):
        pass

    def testContinuousExploration(self):
        etreeSnap = ExplorationTreeSnapshot.load(self.test_files_path + "test-template.xml")
        tree = ExplorationTree.createFromSnapshot(etreeSnap)
        tree.setThreadCount(2)

        for i in range(5):
            tree.generateMorphs()
            tree.sortMorphs()
            tree.filterMorphs(FilterMoprhsOper.ALL)
            tree.extend()
            tree.prune()
            print("Iteration {0}".format(i+1))

        results = 'testing_snapshots'
        run_path_finder(results, self.test_files_path + "test-template.xml", 2)
        shutil.rmtree(results)