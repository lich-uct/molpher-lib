/*
 * File:   MinimalTest.hpp
 * Author: sichom
 *
 * Created on Feb 8, 2016, 11:38:51 AM
 */

#ifndef MINIMALTEST_HPP
#define	MINIMALTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationTree.hpp"
#include "data_structs/ExplorationData.hpp"

#include "SAScore_data_loader.hpp"

#include "selectors/chemoper_selectors.h"
#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"

template<typename Number>
std::string NumberToStr(Number num) {
    std::stringstream ss;
    ss << num;
    return ss.str();
}

void printCandidates(std::shared_ptr<ExplorationTree> tree) {
    int counter = 1;
    for (auto candidate : tree->getCandidateMorphs()) {
        std::cout << NumberToStr(counter++) + ": " << candidate->getSMILES() << " -- " + NumberToStr(candidate->getDistToTarget()) << std::endl;
    }
}

class MinimalTest : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(MinimalTest);

    CPPUNIT_TEST(testMolpherMol);
    CPPUNIT_TEST(testTree);
    CPPUNIT_TEST(testExplorationData);

    CPPUNIT_TEST_SUITE_END();

public:
    MinimalTest();
    virtual ~MinimalTest();
    void setUp();
    void tearDown();

private:
    const std::string test_dir;
    
    void testMolpherMol();
    void testTree();
    void testExplorationData();
};

#endif	/* MINIMALTEST_HPP */

