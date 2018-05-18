"""
This module facilitates access to the random number generator
"""

import molpher.swig_wrappers.core as wrappers

def set_random_seed(seed):
    wrappers.set_random_seed(seed)

def get_random_number(floor, ceil):
    return wrappers.get_random_number(floor, ceil)