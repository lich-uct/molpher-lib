/*
 * File:   MinimalTest.hpp
 * Author: sichom
 *
 * Created on Feb 8, 2016, 11:38:51 AM
 */

#ifndef MINIMALTEST_HPP
#define	MINIMALTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

class MinimalTest : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(MinimalTest);

    CPPUNIT_TEST(testMolpherMol);

    CPPUNIT_TEST_SUITE_END();

public:
    MinimalTest();
    virtual ~MinimalTest();
    void setUp();
    void tearDown();

private:
    void testMolpherMol();
};

#endif	/* MINIMALTEST_HPP */

