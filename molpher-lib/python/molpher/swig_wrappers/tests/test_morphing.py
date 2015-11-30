import unittest

import shutil

from molpher.swig_wrappers.core import *
from pkg_resources import resource_filename

class TestMorphing(unittest.TestCase):

    def setUp(self):
        self.test_files_path = resource_filename('molpher.swig_wrappers.tests', 'test_files/')

    def tearDown(self):
        pass

    def testContinuousExploration(self):
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