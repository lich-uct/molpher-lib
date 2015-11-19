from pkg_resources import Requirement, resource_filename

import molpher.swig_wrappers.core as core

def load_SAScore(path):
    print("Loading data from:", path)
    core.SAScore.loadData(path)

load_SAScore(resource_filename('molpher.swig_wrappers', 'SAScore.dat'))