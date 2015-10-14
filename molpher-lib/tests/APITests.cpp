/*
 * File:   MiscTests.cpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:25 AM
 */

#include "APITests.hpp"

#include "../include/MorphingManager.hpp"
#include "../molpher_API/ExplorationParameters.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(APITests);

APITests::APITests() {
}

APITests::~APITests() {
}

void APITests::setUp() {
}

void APITests::tearDown() {
}

void APITests::testExplorationParametersClass() {
    ExplorationParameters params;
    auto snp = params.createIterationSnapshot();
    CPPUNIT_ASSERT_EQUAL(decltype(snp.jobId)(0), snp.jobId);
    CPPUNIT_ASSERT_EQUAL(decltype(snp.iterIdx)(0), snp.iterIdx);
    CPPUNIT_ASSERT_EQUAL(decltype(snp.elapsedSeconds)(0), snp.elapsedSeconds);
    CPPUNIT_ASSERT(!snp.chemOperSelectors.empty());
}

void APITests::testConstructors() {
    CPPUNIT_ASSERT_THROW(MoprhingManager mng("wrong");, std::runtime_error);
    CPPUNIT_ASSERT_NO_THROW(MoprhingManager mng("../templates/test-template.xml"));
}

