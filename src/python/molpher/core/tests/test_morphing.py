
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

import os, sys
import unittest
from pprint import pprint
from rdkit import Chem

from pkg_resources import resource_filename

from molpher import random
from molpher.examples import morphing_opers, exploration_basics
from molpher.core import ExplorationTree
from molpher.core import MolpherMol
from molpher.core.operations import *

class TestMorphing(unittest.TestCase):

    @staticmethod
    def getPathToMol(tree, mol):
        assert tree.hasMol(mol)
        path = []
        current = mol
        while current.parent_smiles:
            path.append(current)
            current = tree.fetchMol(current.parent_smiles)
        path.reverse()
        return path

    def setUp(self):
        random.set_random_seed(42)
        self.test_source = 'CCO'
        self.test_target = 'C1=COC=C1'
        self.test_dir = os.path.abspath(resource_filename('molpher.core.tests', 'test_files/'))
        self.test_template_path = os.path.join(self.test_dir, 'test-template.xml')
        self.cymene_locked = os.path.join(self.test_dir, 'cymene.sdf')
        self.ethanol_locked = os.path.join(self.test_dir, 'ethanol.sdf')
        self.propanol = os.path.join(self.test_dir, 'propanol.sdf')
        self.remove_bond_test_mol = os.path.join(self.test_dir, 'remove_bond_test_mol.sdf')
        self.alanine = os.path.join(self.test_dir, 'alanine.sdf')
        self.isopropylphenol = os.path.join(self.test_dir, 'isopropylphenol.sdf')
        self.contract_bond_test_mol = os.path.join(self.test_dir, 'contract_bond_test_mol.sdf')
        self.reroute_test_mol = os.path.join(self.test_dir, 'reroute_test_mol.sdf')
        self.captopril = os.path.join(self.test_dir, 'captopril.sdf')

    def tearDown(self):
        pass

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

        my_callback = lambda a, b : a.getDistToTarget() > b.getDistToTarget()
        my_sort = SortMorphsOper(tree, my_callback)
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
        self.assertEqual(len(tree.candidates), selected)

        tree.extend()

        callback = lambda x : sys.stdout.write(x.smiles + ' : ' + str(x.dist_to_target) + '\n')
        oper = TraverseOper(callback=callback)
        tree.runOperation(oper)

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

        all_bad_structures = []
        def collect_nonsyntetizable(morph, operator):
            if morph.sascore > 6:
                all_bad_structures.append(morph)

        class MorphingIteration(TreeOperation):

            parent = self

            def __init__(self, tree):
                super(MorphingIteration, self).__init__()
                self._tree = tree

            def __call__(self):
                print('Iteration: ', self._tree.getGenerationCount() + 1)
                self._tree.generateMorphs([collect_nonsyntetizable])
                for mol in self._tree.candidates:
                    self.parent.assertEqual(None, mol.tree)
                self._tree.sortMorphs()
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
            # , 'threads' : 1
        })

        iterate = MorphingIteration(tree)
        counter = 0
        while True:
            iterate()
            counter += 1
            if tree.path_found:
                target = tree.fetchMol(self.test_target)
                assert target
                print("Path found after {0} iterations:".format(counter))
                path = self.getPathToMol(tree, target)
                pprint([(x.smiles, x.dist_to_target, x.parent_operator) for x in path])
                break

        child = tree.leaves[0]
        self.assertTrue(child.tree)
        self.assertTrue(tree.hasMol(child))
        parent = child.getParentSMILES()
        tree.deleteSubtree(parent)
        self.assertFalse(tree.hasMol(parent))
        self.assertFalse(tree.hasMol(child))
        self.assertEqual(None, child.tree)
        self.assertEqual(parent, child.getParentSMILES())

        # check if valid molecules were extracted
        self.assertTrue(len(all_bad_structures) > 0)
        for mol in all_bad_structures:
            self.assertTrue(mol.smiles)

        # check descendents
        def check_descs(morph):
            for desc_smiles in morph.descendents:
                desc = tree.fetchMol(desc_smiles)
                self.assertTrue(desc.tree)
        tree.traverse(check_descs)

    def testMorphingWithLocks(self):
        tree = ExplorationTree.create(source=MolpherMol(self.captopril))

        # generate two generations of morphs and save them all to a list
        morphs = []
        def some_collector(morph, operator):
            self.assertTrue(operator.name)
            self.assertTrue(morph.smiles)
            morphs.append((morph, operator))
        gen_morphs = GenerateMorphsOper(collectors=[some_collector])
        tree.runOperation(gen_morphs)
        tree.sortMorphs()
        tree.filterMorphs()
        tree.extend()
        tree.runOperation(gen_morphs)
        tree.extend()

        # check if all generated morphs satisfy some conditions
        locked_pattern = Chem.MolFromSmarts('C(=O)N1CCCC1C(=O)O')
        for x in morphs:
            self.assertTrue(x[0].smiles)
            self.assertTrue(x[1].name)
            self.assertTrue(x[0].asRDMol().HasSubstructMatch(locked_pattern))

    def testExamples(self):
        morphing_opers.main(self.captopril)
        exploration_basics.main()

if __name__ == "__main__":
    unittest.main()