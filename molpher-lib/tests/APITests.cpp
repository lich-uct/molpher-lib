/*
 * File:   MiscTests.cpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:25 AM
 */

#include "APITests.hpp"

#include "extensions/SAScore.h"

#include "../include/molpher_API/ExplorationParameters.hpp"
#include "../include/molpher_API/MolpherMol.hpp"
#include "../include/molpher_API/ExplorationTreeSnapshot.hpp"
#include "../include/molpher_API/ExplorationTree.hpp"
#include "molpher_API/operations/FilterMorphsOper.hpp"

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
    CPPUNIT_ASSERT(!params.valid());
    params.setSourceMol("CCO");
    CPPUNIT_ASSERT(params.valid());
}

void APITests::testExplorationTreeSnapshotClass() {
    ExplorationTreeSnapshot* etreeSnap = ExplorationTreeSnapshot::load(test_files_path + "test-template.xml");
    etreeSnap->save(test_files_path + "snappitty_snap.snp");
}

void APITests::testExplorationTreeClass() {
    ExplorationParameters params;
    params.setSourceMol("CCO");
    ExplorationTree param_tree(params);
    ExplorationTreeSnapshot* snap = param_tree.createSnapshot();
    
    std::string smiles("OCCO");
    ExplorationTree smile_tree(smiles);
    snap = smile_tree.createSnapshot();
    snap->save(test_files_path + "snappy_snap.snp");
    snap = snap->load(test_files_path + "snappy_snap.snp");
    ExplorationTree* etree = ExplorationTree::createFromSnapshot(*snap);
    MolpherMol* mol = etree->fetchMol(smiles);
    CPPUNIT_ASSERT_EQUAL(smiles, mol->getSMILES());
}

void APITests::testExploration() {
    SAScore::loadData();
    
    // test morphing methods
    ExplorationTreeSnapshot* etreeSnap = ExplorationTreeSnapshot::load(test_files_path + "test-template.xml");
    ExplorationTree* treept = ExplorationTree::createFromSnapshot(*etreeSnap);
    ExplorationTree& tree = *treept;
    tree.setThreadCount(2);
    
    std::vector<MolpherMol> leaves;
    tree.fetchLeaves(leaves);
    CPPUNIT_ASSERT(1 == leaves.size());
    
    tree.generateMorphs();
    std::vector<MolpherMol> morphs = tree.getCandidateMorphs();
    CPPUNIT_ASSERT(1 == leaves.size());
    CPPUNIT_ASSERT(!morphs.empty());
    for (auto morph : morphs) {
        CPPUNIT_ASSERT(morph.getMol().IsValid());
    }
    
    tree.sortMorphs();
    morphs = tree.getCandidateMorphs();
    MolpherMolecule* previous = nullptr;
    for (auto& morph : morphs) {
        MolpherMolecule& molpher_molecule = morph.getMol();
        CPPUNIT_ASSERT(molpher_molecule.IsValid());
        if (previous) {
            CPPUNIT_ASSERT(molpher_molecule.distToTarget >= previous->distToTarget);
        }
        std::cout << morph.getSMILES() << std::endl;
        previous = &molpher_molecule;
    }
    
    tree.filterMorphs(FilterMoprhsOper::MorphFilters::COUNT | FilterMoprhsOper::MorphFilters::PROBABILITY | FilterMoprhsOper::MorphFilters::WEIGHT);
    std::vector<bool> mask = tree.getCandidateMorphsMask();
    CPPUNIT_ASSERT_EQUAL(mask.size(), morphs.size());
    
    tree.extend(); // TODO: validate that the morphs have really been added to the tree successfully (tested from Python for the time being)
    
    tree.prune(); // TODO: figure out how to test
}

