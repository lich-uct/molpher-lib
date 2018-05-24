"""
Tests registered here should be run when the conda packages are built.
"""

import unittest

TEST_MODULES = [
    # TODO: automatically import all tests from the package
    'molpher.core.tests.test_api',
    'molpher.core.tests.test_morphing',
]

def run(module_list=TEST_MODULES):
    suite = unittest.TestSuite()

    for t in module_list:
        try:
            # If the module defines a suite() function, call it to get the suite.
            mod = __import__(t, globals(), locals(), ['suite'])
            suitefn = getattr(mod, 'suite')
            suite.addTest(suitefn())
        except (ImportError, AttributeError):
            # else, just load all the test cases from the module.
            suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))

    result = unittest.TextTestRunner().run(suite)
    if result.errors:
        exit(1)