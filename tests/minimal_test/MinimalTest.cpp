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
    MolpherMol mol;
    CPPUNIT_ASSERT(true);
}
