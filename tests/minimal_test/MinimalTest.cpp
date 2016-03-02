/*
 * File:   MinimalTest.cpp
 * Author: sichom
 *
 * Created on Feb 8, 2016, 11:38:51 AM
 */
#include <stdexcept>

#include "MinimalTest.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(MinimalTest);

MinimalTest::MinimalTest() :
test_dir("tests/test_files/")
{
    // no action
}

MinimalTest::~MinimalTest() {
    // no action
}

void MinimalTest::setUp() {
    load_data_from("res/SAScore.dat");
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
}

