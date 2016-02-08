/*
 * File:   MinimalTest.cpp
 * Author: sichom
 *
 * Created on Feb 8, 2016, 11:38:51 AM
 */

#include "MinimalTest.hpp"

#include "data_structs/MolpherMol.hpp"


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
    MolpherMol mol("CCO");
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(DBL_MAX, mol.getDistToTarget());
    
    auto mol_copy = mol.copy();
    mol_copy->setDistToTarget(1.0);
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol_copy->getSMILES());
    CPPUNIT_ASSERT_EQUAL(1.0, mol_copy->getDistToTarget());
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(DBL_MAX, mol.getDistToTarget());
    
    CPPUNIT_ASSERT(!mol.getTree());
    CPPUNIT_ASSERT(!mol_copy->getTree());
}
