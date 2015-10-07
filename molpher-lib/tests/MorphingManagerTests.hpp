/*
 * File:   MiscTests.hpp
 * Author: sichom
 *
 * Created on Oct 7, 2015, 9:45:24 AM
 */

#ifndef MISCTESTS_HPP
#define	MISCTESTS_HPP

#include <cppunit/extensions/HelperMacros.h>

class MorphingManagerTests : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(MorphingManagerTests);

    CPPUNIT_TEST(testConstructors);

    CPPUNIT_TEST_SUITE_END();

public:
    MorphingManagerTests();
    virtual ~MorphingManagerTests();
    void setUp();
    void tearDown();

private:
    void testConstructors();
};

#endif	/* MISCTESTS_HPP */

