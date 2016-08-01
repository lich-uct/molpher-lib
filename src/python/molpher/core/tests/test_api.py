
# Copyright (c) 2016 Martin Sicho
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import unittest

from pkg_resources import resource_filename

from molpher import random
from molpher.core import ExplorationTree
from molpher.core import MolpherMol
from molpher.core.operations import *
from molpher.core.operations.callbacks import SortMorphsCallback
from molpher.core.selectors import *
from molpher.core import ExplorationData

random.set_random_seed(42)

class TestPythonAPI(unittest.TestCase):
    
    def setUp(self):
        self.test_source = 'CCO'
        self.test_target = 'C1=COC=C1'
        self.test_dir = os.path.abspath(resource_filename('molpher.core.tests', 'test_files/'))
        self.test_template_path = os.path.join(self.test_dir, 'test-template.xml')

    def tearDown(self):
        pass

    def testMolpherMol(self):
        mol = MolpherMol(self.test_target)
        mol.smiles = 'CCC'
        self.assertTrue('CCC')
        # mol.historic_descendents = ('CCC', 'CCCC')
        # self.assertEqual(('CCC', 'CCCC'), mol.historic_descendents)

        copy = mol.copy()
        copy.sascore = 0.54
        self.assertEqual(0.54, copy.sascore)
        self.assertEqual(0, mol.sascore)

        tree = ExplorationTree.create(source=mol.smiles, target='CCCNCCC')
        tree = ExplorationTree.create(source=mol, target='CCCNCCC')
        tree = ExplorationTree.create(source=mol, target=MolpherMol('CCCNCCC'))
        self.assertTrue(tree.hasMol(mol))
        def assign(x):
            tree.fetchMol(mol.smiles).smiles = x
        self.assertRaises(RuntimeError, assign, 'CCO')

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
        tree_from_dict = ExplorationTree.create(tree_data=params_dict)
        tree_from_params = ExplorationTree.create(tree_data=params)
        tree_from_SMILES = ExplorationTree.create(source=mol1, target=mol2)
        def test_tree(tree):
            self.assertEqual(tree.params['source'], mol1)
            self.assertEqual(tree.params['target'], mol2)

        test_tree(tree_from_dict)
        test_tree(tree_from_params)
        test_tree(tree_from_SMILES)

        tree = tree_from_params

        # if we try to set source for non-empty tree, exception should be raised
        def func():
            tree.params = {
                'source' : mol2
                , 'target' : 'C'
            }
        self.assertRaises(RuntimeError, func)

        tree.thread_count = 1
        tree.params = {
            'target' : 'C'
        }
        self.assertEqual(1, tree.thread_count)
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], 'C')
        self.assertEqual(tree.params['operators'], params.param_dict['operators']) # we should still have the same opers set

        tree.params = params; tree.thread_count = 0 # assign the original parameters back
        self.assertEqual(0, tree.thread_count)
        self.assertEqual(tree.params['source'], mol1)
        self.assertEqual(tree.params['target'], mol2)
        self.assertEqual(tree.params['operators'], params.param_dict['operators'])

        leaf = tree.leaves[0]
        self.assertRaises(RuntimeError, lambda : leaf.setSMILES('CCCC'))
        self.assertTrue(tree.hasMol(leaf))
        # self.assertEqual(tree, leaf.tree) # FIXME: add a reliable operator for comparison between trees
        leaf.setDistToTarget(0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)

        leaf_copy = tree.leaves[0].copy()
        # self.assertFalse(tree.hasMol(leaf_copy)) # FIXME: add a reliable operator for comparison between trees (this should check both the SMILES and the tree ownership)
        self.assertEqual(leaf_copy.getDistToTarget(), 0.5)
        leaf_copy.setDistToTarget(0.7)
        self.assertEqual(leaf.getDistToTarget(), 0.5)
        self.assertEqual(tree.leaves[0].getDistToTarget(), 0.5)
        self.assertEqual(leaf_copy.getDistToTarget(), 0.7)

    def testOperations(self):
        tree = ExplorationTree.create(tree_data={
            'source' : self.test_source
            , 'target' : self.test_target
        })

        iteration = [
            GenerateMorphsOper()
            , SortMorphsOper()
            , FilterMorphsOper()
            , ExtendTreeOper()
            , PruneTreeOper()
        ]

        for oper in iteration:
            self.assertRaises(RuntimeError, lambda : oper())

        fl = FindLeavesOper()
        for oper in iteration:
            tree.runOperation(oper)
        tree.runOperation(fl)
        for leaf1, leaf2, leaf3 in zip(sorted(fl.leaves), sorted(fl.tree.leaves), sorted(tree.leaves)):
            self.assertTrue(leaf1.smiles == leaf2.smiles == leaf3.smiles)

        tree.generateMorphs()
        tree.sortMorphs()
        previous = None
        for morph in tree.candidates:
            if previous:
                self.assertTrue(morph.dist_to_target >= previous)
                previous = morph.dist_to_target
            else:
                previous = morph.dist_to_target
        print([x.dist_to_target for x in tree.candidates])

        class MySort(SortMorphsCallback):

            def __call__(self, a, b):
                return a.getDistToTarget() > b.getDistToTarget()

        my_callback = MySort()
        my_sort = SortMorphsOper(tree, my_callback) # x(tree, MySort()) gives a segfault
        my_sort()

        previous = None
        for morph in tree.candidates:
            if previous:
                self.assertTrue(morph.dist_to_target <= previous)
                previous = morph.dist_to_target
            else:
                previous = morph.dist_to_target
        print([x.dist_to_target for x in tree.candidates])

        tree.filterMorphs()
        selected = sum(tree.candidates_mask)
        clean_stuff = CleanMorphsOper()
        tree.runOperation(clean_stuff)
        self.assertEquals(len(tree.candidates), selected)

    def testMorphing(self):
        def callback(morph):
            callback.morphs_in_tree += 1
            self.assertTrue(morph)
            self.assertTrue(morph.tree)
            if morph.getItersWithoutDistImprovement() > 3:
                print('Callback output:')
                print(morph.getSMILES(), morph.getItersWithoutDistImprovement(), morph.getDistToTarget())
            if not callback.closest_mol:
                callback.closest_mol = morph
            current_dist = morph.getDistToTarget()
            min_dist = callback.closest_mol.getDistToTarget()
            if min_dist > current_dist:
                callback.closest_mol = morph
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

        tree = ExplorationTree.create(tree_data={
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

        # check descendents
        def check_descs(morph):
            for desc_smiles in morph.descendents:
                desc = tree.fetchMol(desc_smiles)
                self.assertTrue(desc.tree)
        tree.traverse(check_descs)

if __name__ == "__main__":
    unittest.main()