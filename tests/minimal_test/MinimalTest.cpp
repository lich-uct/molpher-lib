/*
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdexcept>
#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

#include "MinimalTest.hpp"
#include "io/stdout.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(MinimalTest);

MinimalTest::MinimalTest() :
test_dir("../test_files/")
{
    // no action
}

MinimalTest::~MinimalTest() {
    // no action
}

void MinimalTest::setUp() {
    load_data_from("../../res/SAScore.dat");
}

void MinimalTest::tearDown() {
    // no action
}

void MinimalTest::testMolpherMol() {
    // test the simple constructor
    MolpherMol mol("CCO");
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(DBL_MAX, mol.getDistToTarget());
    
    // test the copy method
    auto mol_copy = mol.copy();
    mol_copy->setDistToTarget(1.0);
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol_copy->getSMILES());
    CPPUNIT_ASSERT_EQUAL(1.0, mol_copy->getDistToTarget());
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(DBL_MAX, mol.getDistToTarget());
    
    // molecules created by themselves should not belong to a tree
    CPPUNIT_ASSERT(!mol.getTree());
    CPPUNIT_ASSERT(!mol_copy->getTree());
    CPPUNIT_ASSERT(!mol.isBoundToTree());
    CPPUNIT_ASSERT(!mol_copy->isBoundToTree());
    
    mol.setSMILES("C1CCC1");
    CPPUNIT_ASSERT_EQUAL(std::string("C1CCC1"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol_copy->getSMILES());
    
//    CPPUNIT_ASSERT_THROW(mol.setSMILES("ABCDE");, RDKit::SmilesParseException); // rdkit also gives a segmentation fault for some reason...
    
    // test the empty constructor (molecule with empty SMILES should not be valid and all values should be initialized to defaults)
    MolpherMol empty;
    CPPUNIT_ASSERT_EQUAL(std::string(""), empty.getSMILES());
    CPPUNIT_ASSERT(!empty.isValid());
    CPPUNIT_ASSERT_EQUAL(OP_ADD_ATOM, static_cast<ChemOperSelector>(empty.getParentOper()));
    CPPUNIT_ASSERT_EQUAL((unsigned) 0, empty.getItersWithoutDistImprovement());
    
    // test the complete constructor
    MolpherMol complete("CC(=O)C", "C3O", "NC(=O)C",
                OP_MUTATE_ATOM, 0.5, 0.0,
                0.0, 3.1);
    CPPUNIT_ASSERT_EQUAL(OP_MUTATE_ATOM, static_cast<ChemOperSelector>(complete.getParentOper()));
    CPPUNIT_ASSERT_EQUAL(std::string("NC(=O)C"), complete.getParentSMILES());
}

void MinimalTest::testExplorationData() {
    // test empty initialzation, default data should not be valid
    ExplorationData data;
    CPPUNIT_ASSERT(!data.isValid());
    
    // get the source molecule, should have no SMILES
    auto source = data.getSource();
    CPPUNIT_ASSERT_EQUAL(std::string(""), source->getSMILES());
    
    // modify the source and target molecule and save them back to tree data
    source->setSMILES("NC(=O)C");
    auto target = data.getTarget();
    target->setSMILES("CC(=O)C");
    data.setSource(*source);
    data.setTarget(*target);
    
    // data should now be valid
    CPPUNIT_ASSERT(data.isValid());
    
    // the tree should now contain the source as root
    auto map = data.getTreeMap();
    CPPUNIT_ASSERT(map->size() == 1);
    CPPUNIT_ASSERT_EQUAL(map->find(source->getSMILES())->second->getSMILES(), source->getSMILES());
    
    // we shouldn't be able to set source molecule for the second time
    CPPUNIT_ASSERT_THROW(data.setSource(*source);, std::runtime_error);
    data.setTarget(*target); // we can change the target though
    
    // check the selectors
    CPPUNIT_ASSERT_EQUAL(FP_MORGAN, static_cast<FingerprintSelector>(data.getFingerprint()));
    CPPUNIT_ASSERT_EQUAL(SC_TANIMOTO, static_cast<SimCoeffSelector>(data.getSimilarityCoefficient()));
    auto operators = data.getChemicalOperators();
    CPPUNIT_ASSERT(!operators.empty());
    for (auto& oper : operators) {
        std::cout << ChemOperLongDesc(oper) << std::endl;
    }
    
    // test the serialization
    auto data_from_template = ExplorationData::load(test_dir + "test-template.xml");
    CPPUNIT_ASSERT(data_from_template->isValid());
    
    data_from_template->setTarget(*target);
    data_from_template->save(test_dir + "template_snapshot.xml");
    auto data_from_snapshot = ExplorationData::load(test_dir + "template_snapshot.xml");
    CPPUNIT_ASSERT_EQUAL(target->getSMILES(), data_from_snapshot->getTarget()->getSMILES());
}


void MinimalTest::testTree() {
    // tree creation
    auto tree = ExplorationTree::create("CCO", "C1CCC1");
    CPPUNIT_ASSERT_THROW(ExplorationTree::create("CCO", "");, std::runtime_error);
    
    // retrieve the source and see if it belongs to the correct tree
    auto source = tree->fetchMol("CCO");
    CPPUNIT_ASSERT(source->isBoundToTree());
    auto tree_of_source = source->getTree();
    CPPUNIT_ASSERT(tree_of_source);
    CPPUNIT_ASSERT_EQUAL(tree, tree_of_source);
    CPPUNIT_ASSERT(tree_of_source->hasMol(source));
    CPPUNIT_ASSERT(tree_of_source->hasMol(source->getSMILES()));
    
    // we should not be able to set ownership of an already owned molecule
    CPPUNIT_ASSERT_THROW(source->setOwner(tree);, std::runtime_error);
    
    // fetch leaves shouldn't increase the distance imporvement counter for source
    auto leaves = tree->fetchLeaves(true);
    CPPUNIT_ASSERT(leaves.size() == 1);
    CPPUNIT_ASSERT_EQUAL(source, leaves[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned) 0, source->getItersWithoutDistImprovement());
    
    // generate some candidate morphs
    CPPUNIT_ASSERT(tree->getCandidateMorphs().empty());
    CPPUNIT_ASSERT(tree->getCandidateMorphsMask().empty());
    tree->generateMorphs();
    CPPUNIT_ASSERT(!tree->getCandidateMorphs().empty());
    CPPUNIT_ASSERT(!tree->getCandidateMorphsMask().empty());
    CPPUNIT_ASSERT_EQUAL(tree->getCandidateMorphs().size(), tree->getCandidateMorphsMask().size());
    for (auto candidate : tree->getCandidateMorphs()) {
        CPPUNIT_ASSERT(!candidate->isBoundToTree());
        CPPUNIT_ASSERT_EQUAL(false, (bool) candidate->getTree());
    }
    printCandidates(tree);
    
    // sort the morphs
    tree->sortMorphs();
    printCandidates(tree);
    double previous = 0;
    for (auto candidate : tree->getCandidateMorphs()) {
        CPPUNIT_ASSERT(candidate->getDistToTarget() >= previous);
        previous = candidate->getDistToTarget();
    }
    
    // filter morphs
    tree->filterMorphs(true);
    printCandidates(tree);
    
    // extend the tree
    CPPUNIT_ASSERT_EQUAL((unsigned) 0, tree->getGenerationCount());
    
    tree->extend();
    leaves = tree->fetchLeaves();
    auto source_descendents = source->getDescendants();
    for (auto leaf : leaves) {
        CPPUNIT_ASSERT(leaf->isBoundToTree());
        CPPUNIT_ASSERT_EQUAL(tree, leaf->getTree());
        
        auto parent_smiles = leaf->getParentSMILES();
        CPPUNIT_ASSERT_EQUAL(source->getSMILES(), parent_smiles);
        CPPUNIT_ASSERT_EQUAL(source, tree->fetchMol(parent_smiles));
        
        auto it = source_descendents.find(leaf->getSMILES());
        CPPUNIT_ASSERT(it != source_descendents.end());
        CPPUNIT_ASSERT_EQUAL(leaf, tree->fetchMol(*it));
    }
    
    CPPUNIT_ASSERT(tree->getCandidateMorphs().empty());
    CPPUNIT_ASSERT(tree->getCandidateMorphsMask().empty());
    CPPUNIT_ASSERT_EQUAL((unsigned) 1, tree->getGenerationCount());
    
    // prune the tree
    tree->prune();
    auto leaf_zero = leaves[0];
    tree->deleteSubtree(leaf_zero->getSMILES());
    CPPUNIT_ASSERT(!tree->hasMol(leaf_zero->getSMILES()));
    CPPUNIT_ASSERT(!leaf_zero->getTree());
    leaf_zero->removeFromTree(); // shouldn't do anything
    
    // remove all newly generated morphs
    CPPUNIT_ASSERT_THROW(tree->deleteSubtree(source->getSMILES(), false);, std::runtime_error); // should not be possible
    tree->deleteSubtree(source->getSMILES(), true); // remove only descendents
    
    source_descendents = source->getDescendants();
    auto source_hist_descendents = source->getHistoricDescendants();
    for (auto leaf : leaves) {
        CPPUNIT_ASSERT(!tree->hasMol(leaf->getSMILES()));
        CPPUNIT_ASSERT(!leaf->getTree());
        
        auto it = source_descendents.find(leaf->getSMILES());
        CPPUNIT_ASSERT(it == source_descendents.end());
        CPPUNIT_ASSERT(it != source_hist_descendents.end());
    }
    
    // make a few more generations
    for (unsigned iter_idx = 0; iter_idx != 5; iter_idx++) {
        tree->generateMorphs();
        tree->sortMorphs();
        tree->filterMorphs();
        printCandidates(tree);
        tree->extend();
        std::cout << "Path found: " + NumberToStr(tree->isPathFound()) << std::endl;
        tree->prune();
    }
    
    // test the tree traversal
    PrintMols printing_callback;
    tree->traverse(printing_callback);
    
    // serialize the tree into file
    tree->save(test_dir + "testTree_snapshot.xml");
    
    // create new tree from the file
    auto tree_from_file = ExplorationTree::create(test_dir + "testTree_snapshot.xml");
    
    // make a few more generations on the loaded tree
    for (unsigned iter_idx = 0; iter_idx != 3; iter_idx++) {
        tree_from_file->generateMorphs();
        tree_from_file->sortMorphs();
        tree_from_file->filterMorphs();
        printCandidates(tree_from_file);
        tree_from_file->extend();
        std::cout << "Path found: " + NumberToStr(tree_from_file->isPathFound()) << std::endl;
        tree_from_file->prune();
    }
}

void MinimalTest::testRDKit() {
    //TODO: remove this from tests when the fixed_atom feature is ready

    RDKit::ROMol* mol = RDKit::SDMolSupplier(test_dir + "Structure2D_CID_4914.sdf").next();
    print_mol_info(mol);

//    RDKit::RWMol* morphed_mol = RDKit::SmilesToMol(RDKit::MolToSmiles(*mol));
    RDKit::RWMol morphed_mol(*mol);
    RDKit::MolOps::Kekulize(morphed_mol);
    print_mol_info((RDKit::ROMol*) &morphed_mol);

//    morphed_mol.removeAtom((uint) 0);
//    print_mol_info((RDKit::ROMol*) &morphed_mol);

    uint bindingAtomIdx = 3;
    RDKit::Atom* bindingAtom = morphed_mol.getAtomWithIdx(bindingAtomIdx);
    RDKit::Atom addedAtom = *(morphed_mol.getAtomWithIdx(bindingAtomIdx));
    uint newAtomIdx = morphed_mol.addAtom(&addedAtom);
    morphed_mol.addBond(bindingAtom->getIdx(), newAtomIdx, RDKit::Bond::SINGLE);
    print_mol_info((RDKit::ROMol*) &morphed_mol);
}

