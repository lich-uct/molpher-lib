/*
 * File:   MiscTests.hpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:24 AM
 */

#ifndef MISCTESTS_HPP
#define	MISCTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <chrono>
#include <thread>

#include "molpher_API/callbacks/TraverseCallback.hpp"

class PrintMolCallback : public TraverseCallback {
    
    public:
        PrintMolCallback() : TraverseCallback() {
            // no action
        }

        virtual void processMorph(MolpherMol& morph) {
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
            std::cout << morph.getSMILES() + " -- parent: " + morph.fetchMolpherMolecule().parentSmile << std::endl;
        }

};

class APITests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(APITests);

    CPPUNIT_TEST(testExplorationParametersClass);
    CPPUNIT_TEST(testMolpherMolClass);
    CPPUNIT_TEST(testExplorationTreeSnapshotClass);
    CPPUNIT_TEST(testExplorationTreeClass);
    CPPUNIT_TEST(testExploration);

    CPPUNIT_TEST_SUITE_END();

public:
    APITests();
    virtual ~APITests();
    void setUp();
    void tearDown();

private:
    std::string test_files_path;
    
    void testMolpherMolClass();
    void testExplorationParametersClass();
    void testExplorationTreeSnapshotClass();
    void testExplorationTreeClass();
    void testExploration();
};

#endif	/* MISCTESTS_HPP */

