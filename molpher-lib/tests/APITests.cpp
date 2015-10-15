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
#include "../include/molpher_API/ExplorationTree.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(APITests);

APITests::APITests() : test_files_path("tests/test_files/") {
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
    params.setSourceMol("CCO");
    CPPUNIT_ASSERT(params.valid());
}

void APITests::testExplorationTreeSnapshotClass() {
    ExplorationTreeSnapshot etreeSnap = ExplorationTreeSnapshot::load(test_files_path + "test-template.xml");
    etreeSnap.save(test_files_path + "snappitty_snap.snp");
}

void APITests::testExplorationTreeClass() {
    ExplorationParameters params;
    params.setSourceMol("CCO");
    ExplorationTree param_tree(params);
    ExplorationTree smile_tree("OCCO");
    ExplorationTreeSnapshot snap = smile_tree.createSnapshot();
    snap.save(test_files_path + "snappy_snap.snp");
    snap = snap.load(test_files_path + "snappy_snap.snp");
    ExplorationTree file_tree = ExplorationTree::createFromSnapshot(snap);
    file_tree.setThreadCount(2);
    std::vector<MolpherMol> ret = file_tree.fetchLeaves();
    CPPUNIT_ASSERT(!ret.empty());
}

