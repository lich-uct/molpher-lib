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

#ifndef MINIMALTEST_HPP
#define	MINIMALTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationTree.hpp"
#include "data_structs/ExplorationData.hpp"
#include "operations/callbacks/TraverseCallback.hpp"

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
    auto mask = tree->getCandidateMorphsMask();
    for (auto candidate : tree->getCandidateMorphs()) {
        bool mask_val = mask[counter - 1];
        std::cout 
                << NumberToStr(counter++) + ": " 
                << candidate->getSMILES() << " -- " 
                << NumberToStr(candidate->getDistToTarget()) 
                << "(" + NumberToStr(mask_val) + ")"
                << std::endl;
    }
}

class PrintMols : public TraverseCallback {
    
    virtual void operator()(std::shared_ptr<MolpherMol> morph) const {
        std::cout << morph->getSMILES() << " -- descendants:" << std::endl;
        for (auto descendant : morph->getDescendants()) {
            std::cout << "\t" << descendant << std::endl;
        }
    }
    
public:
    
    PrintMols() : TraverseCallback() {
        // no action
    }

};

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

