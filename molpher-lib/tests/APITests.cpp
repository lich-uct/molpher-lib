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
#include "molpher_API/operations/TraverseOper.hpp"

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
    MolpherMol mol;
    CPPUNIT_ASSERT_EQUAL(decltype(mol.fetchMolpherMolecule().itersWithoutDistImprovement)(0), mol.fetchMolpherMolecule().itersWithoutDistImprovement);
    CPPUNIT_ASSERT_EQUAL(decltype(mol.fetchMolpherMolecule().smile)(""), mol.fetchMolpherMolecule().smile);
    MolpherMol mol2;
    mol2.fetchMolpherMolecule().smile = "CC";
    mol = mol2;
    CPPUNIT_ASSERT_EQUAL(mol.getSMILES(), mol2.getSMILES());
    
    MolpherMol mol3(mol2);
    CPPUNIT_ASSERT_EQUAL(mol2.getSMILES(), mol3.getSMILES());
    
    CPPUNIT_ASSERT((mol.isBound() && mol2.isBound() && mol3.isBound()) == false);
    
    MolpherMol* mol4 = mol.copy();
    CPPUNIT_ASSERT_EQUAL(mol.getSMILES(), mol4->getSMILES());
    CPPUNIT_ASSERT(&mol != mol4);
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
    ExplorationParameters& generated_params = param_tree.getParams();
    MolpherMol* source = generated_params.getSourceMol();
    MolpherMol* target = generated_params.getTargetMol();
    delete &generated_params;
    std::string ssmiles = source->getSMILES();
    std::string tsmiles = target->getSMILES();
    
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
        CPPUNIT_ASSERT(morph.fetchMolpherMolecule().IsValid());
    }
    
    tree.sortMorphs();
    morphs = tree.getCandidateMorphs();
    MolpherMolecule* previous = nullptr;
    for (auto& morph : morphs) {
        MolpherMolecule& molpher_molecule = morph.fetchMolpherMolecule();
        CPPUNIT_ASSERT(molpher_molecule.IsValid());
        if (previous) {
            CPPUNIT_ASSERT(molpher_molecule.distToTarget >= previous->distToTarget);
        }
        std::cout << morph.getSMILES() << std::endl;
        previous = &molpher_molecule;
    }
    
    tree.filterMorphs(FilterMorphsOper::MorphFilters::ALL);
    std::vector<bool> mask = tree.getCandidateMorphsMask();
    CPPUNIT_ASSERT_EQUAL(mask.size(), morphs.size());
    
    tree.extend();
    
    PrintMolCallback callback;
    TraverseOper traverse(tree, callback);
    int counter = 0;
    while (counter < 20) {
        traverse();
        counter++;
    }
    
    tree.prune();
}

void APITests::testMisc() {
    ExplorationParameters params;
    params.setSourceMol("CCO");
    params.setTargetMol("O1C=CC=C1");
    ExplorationTree tree(params);
    
    const std::vector<MolpherMol>& leaves = tree.fetchLeaves();
    CPPUNIT_ASSERT(leaves.size() > 0);
    for (auto& mol : leaves) {
        CPPUNIT_ASSERT(mol.isBound());
    }
    
    tree.generateMorphs();
    const std::vector<MolpherMol>& candidates = tree.getCandidateMorphs();
    CPPUNIT_ASSERT(candidates.size() > 0);
    for (auto& mol : candidates) {
        CPPUNIT_ASSERT(mol.isBound());
    }
    tree.extend();
    
    delete &leaves;
    leaves = tree.fetchLeaves();
    CPPUNIT_ASSERT(leaves.size() > 0);
    for (auto& mol : leaves) {
        CPPUNIT_ASSERT(mol.isBound());
    }
}

