/*
 * File:   MiscTests.cpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:25 AM
 */

#include "APITests.hpp"

#include "../include/molpher_API/ExplorationParameters.hpp"
#include "../include/molpher_API/MolpherMol.hpp"
#include "../include/molpher_API/ExplorationTreeSnapshot.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(APITests);

APITests::APITests() {
}

APITests::~APITests() {
}

void APITests::setUp() {
}

void APITests::tearDown() {
}

void APITests::testMolpherMolClass() {
    MolpherMol mol("CCO");
    CPPUNIT_ASSERT_EQUAL(decltype(mol.getMol().itersWithoutDistImprovement)(0), mol.getMol().itersWithoutDistImprovement);
    CPPUNIT_ASSERT_EQUAL(decltype(mol.getMol().smile)("CCO"), mol.getMol().smile);
}

void APITests::testExplorationParametersClass() {
    ExplorationParameters params;
    auto snp = params.createIterationSnapshot();
    CPPUNIT_ASSERT_EQUAL(decltype(snp.jobId)(0), snp.jobId);
    CPPUNIT_ASSERT_EQUAL(decltype(snp.iterIdx)(0), snp.iterIdx);
    CPPUNIT_ASSERT_EQUAL(decltype(snp.elapsedSeconds)(0), snp.elapsedSeconds);
    CPPUNIT_ASSERT(!snp.chemOperSelectors.empty());
}

void APITests::testExplorationTreeSnapshot() {
    IterationSnapshot snp;
    CPPUNIT_ASSERT_THROW(ExplorationTreeSnapshot esnp(snp), std::runtime_error);
    ExplorationTreeSnapshot etreeSnap = ExplorationTreeSnapshot::load("../templates/test-template.xml");
}

