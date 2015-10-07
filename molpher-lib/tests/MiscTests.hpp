/*
 * File:   MiscTests.hpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:24 AM
 */

#ifndef MISCTESTS_HPP
#define	MISCTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

class MiscTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(MiscTests);

    CPPUNIT_TEST(testMethod);
    CPPUNIT_TEST(testFailedMethod);

    CPPUNIT_TEST_SUITE_END();

public:
    MiscTests();
    virtual ~MiscTests();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();
};

#endif	/* MISCTESTS_HPP */

