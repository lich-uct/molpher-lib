/*
 * File:   MiscTests.hpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:24 AM
 */

#ifndef MISCTESTS_HPP
#define	MISCTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

class APITests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(APITests);

    CPPUNIT_TEST(testConstructors);
    CPPUNIT_TEST(testExplorationParametersClass);
    CPPUNIT_TEST(testMolpherMolClass);
    CPPUNIT_TEST(testExplorationTreeSnapshot);

    CPPUNIT_TEST_SUITE_END();

public:
    APITests();
    virtual ~APITests();
    void setUp();
    void tearDown();

private:
    void testMolpherMolClass();
    void testExplorationParametersClass();
    void testConstructors();
    void testExplorationTreeSnapshot();
};

#endif	/* MISCTESTS_HPP */

