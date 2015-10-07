/*
 * File:   MiscTests.cpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:25 AM
 */

#include "MorphingManagerTests.hpp"

#include "../include/MorphingManager.hpp"


CPPUNIT_TEST_SUITE_REGISTRATION(MorphingManagerTests);

MorphingManagerTests::MorphingManagerTests() {
}

MorphingManagerTests::~MorphingManagerTests() {
}

void MorphingManagerTests::setUp() {
}

void MorphingManagerTests::tearDown() {
}

void MorphingManagerTests::testConstructors() {
    CPPUNIT_ASSERT_THROW(MoprhingManager mng("wrong");, std::runtime_error);
    CPPUNIT_ASSERT_NO_THROW(MoprhingManager mng("../templates/test-template.xml"));
}

