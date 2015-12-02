from pkg_resources import resource_filename

import molpher.swig_wrappers.core as wrappers
from molpher import core

def load_SAScore(path):
    print("Loading data from:", path)
    wrappers.SAScore.loadData(path)

load_SAScore(resource_filename('molpher.swig_wrappers', 'SAScore.dat'))