import os
import sys
from logging import warning

from molpher import random

from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

# set seed to a fixed value (None or False for a random seed)
SEED = None

# dir for stored data
STORAGE_DIR = os.path.abspath('data')

# number of threads to use
MAX_THREADS = 4

# setup pharmacophore fingerprints
FG_BINS = [(0, 2), (2, 5), (5, 8)] # distance bins
FG_MIN_POINTS = 2 # min number of features
FG_MAX_POINTS = 3 # max number of features

# threshold for the bits in common percentage
COMMON_BITS_PERC_THRS = 0.8

# number of iterations to wait before the antidecoys filter is applied
WAIT_FOR_ANTIDECOYS = 3

# turn antidecoys filter off when the distance between two closest molecules between the two trees gets below this value
ANTIDECOYS_DISTANCE_SWITCH = 0.2

# maximum number of paths to find
PATHS_TO_FIND = 100

# maximum iterations to spend on one path in seconds
MAX_ITERS_PER_PATH = 100

def init():
    thismodule = sys.modules[__name__]
    if not hasattr(thismodule, 'INITIALIZED'):
        if SEED not in (None, False):
            random.set_random_seed(SEED)

        if not os.path.exists(STORAGE_DIR):
            os.mkdir(STORAGE_DIR)

        # important stuff for the pharmacophore fingerprints
        setattr(thismodule, 'FDEF_FILE', os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')) # get basic feature definitions
        setattr(thismodule, 'FEATURE_FACTORY', ChemicalFeatures.BuildFeatureFactory(FDEF_FILE)) # make feature factory
        sig_fac = SigFactory(FEATURE_FACTORY, minPointCount=FG_MIN_POINTS, maxPointCount=FG_MAX_POINTS, trianglePruneBins=False)  # make signature factory
        sig_fac.SetBins(FG_BINS)
        sig_fac.Init()
        setattr(thismodule, 'SIG_FAC', sig_fac)

        setattr(thismodule, 'INITIALIZED', True)
    elif not thismodule.INITIALIZED:
        delattr(thismodule, 'INITIALIZED')
        init()
    else:
        warning('Settings module already initialized. Skipping...')

if not hasattr(sys.modules[__name__], 'INITIALIZED'):
    init()