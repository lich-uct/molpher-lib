/*
 * File:   MinimalTest.cpp
 * Author: sichom
 *
 * Created on Feb 8, 2016, 11:38:51 AM
 */
#include <stdexcept>

#include <GraphMol/SmilesParse/SmilesParse.h>

#include "MinimalTest.hpp"

#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationTree.hpp"
#include "data_structs/ExplorationData.hpp"

#include "selectors/chemoper_selectors.h"
#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"


CPPUNIT_TEST_SUITE_REGISTRATION(MinimalTest);

MinimalTest::MinimalTest() {
}

MinimalTest::~MinimalTest() {
}

void MinimalTest::setUp() {
}

void MinimalTest::tearDown() {
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
}


void MinimalTest::testTree() {
    ExplorationTree tree("CCO", "C1CCC1");
    CPPUNIT_ASSERT_THROW(ExplorationTree tree_no_target("CCO", "");, std::runtime_error);
}

