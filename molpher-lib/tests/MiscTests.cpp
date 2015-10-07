/*
 * File:   MiscTests.cpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:25 AM
 */

#include "MiscTests.hpp"


CPPUNIT_TEST_SUITE_REGISTRATION(MiscTests);

MiscTests::MiscTests() {
}

MiscTests::~MiscTests() {
}

void MiscTests::setUp() {
}

void MiscTests::tearDown() {
}

void MiscTests::testMethod() {
    CPPUNIT_ASSERT(true);
}

void MiscTests::testFailedMethod() {
    CPPUNIT_FAIL("this test always fails...");
}

