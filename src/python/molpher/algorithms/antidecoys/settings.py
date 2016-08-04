def init():
    import sys

    # fetch this module
    THIS_MODULE = sys.modules[__name__]

    if not hasattr(THIS_MODULE, 'INITIALIZED'):
        import os

        from molpher import random

        from rdkit import RDConfig
        from rdkit.Chem import ChemicalFeatures
        from rdkit.Chem.Pharm2D.SigFactory import SigFactory

        # set seed to a fixed value (None or False for a random seed)
        THIS_MODULE.SEED = None

        # dir for stored data
        THIS_MODULE.STORAGE_DIR = os.path.abspath('data')

        # number of threads to use
        THIS_MODULE.MAX_THREADS = 4

        # setup pharmacophore fingerprints
        THIS_MODULE.FG_BINS = [(0, 2), (2, 5), (5, 8)] # distance bins
        THIS_MODULE.FG_MIN_POINTS = 2 # min number of features
        THIS_MODULE.FG_MAX_POINTS = 3 # max number of features

        # top farthest molecules from antidecoys to select
        THIS_MODULE.MAX_ANTIFP_SURVIVORS = 500

        # threshold for the mean bits in common percentage
        THIS_MODULE.COMMON_BITS_PERC_THRS = 0.75

        # maximum acceptable common bits percentage in accapted morphs
        THIS_MODULE.COMMON_BITS_PERC_THRS_MIN = 0.75

        # minimum number of iterations where antidecoys are optimized
        THIS_MODULE.MIN_ANTIDECOY_ITERS = 10
        THIS_MODULE.MAX_ANTIDECOY_ITERS = 50

        # paths containing molecules within this distance from target will be rolled back ROLLBACK_MAX_ITERS_ON_CLOSEST iterations during reset
        THIS_MODULE.RESET_CLOSEST_THRESHOLD = 0.5

        # maximum number of iterations to roll back an unsuitable path
        THIS_MODULE.ROLLBACK_MAX_ITERS = 5

        # maximum number of iterations to roll back when removing the closest molecules
        THIS_MODULE.ROLLBACK_MAX_ITERS_ON_CLOSEST = 1

        # number of iterations to wait before the antidecoys filter is applied
        THIS_MODULE.WAIT_FOR_ANTIDECOYS = 0

        # turn antidecoys filter off when the distance between two closest molecules between the two trees gets below this value
        THIS_MODULE.ANTIDECOYS_DISTANCE_SWITCH = 0.4

        # maximum number of paths to find
        THIS_MODULE.PATHS_TO_FIND = 50

        # maximum iterations to spend on one path in seconds
        THIS_MODULE.MAX_ITERS_PER_PATH = 100

        if THIS_MODULE.SEED not in (None, False):
            random.set_random_seed(THIS_MODULE.SEED)

        if not os.path.exists(THIS_MODULE.STORAGE_DIR):
            os.mkdir(THIS_MODULE.STORAGE_DIR)

        # important stuff for the pharmacophore fingerprints
        THIS_MODULE.FDEF_FILE = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef') # get basic feature definitions
        THIS_MODULE.FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(THIS_MODULE.FDEF_FILE) # make feature factory
        sig_fac = SigFactory(THIS_MODULE.FEATURE_FACTORY, minPointCount=THIS_MODULE.FG_MIN_POINTS, maxPointCount=THIS_MODULE.FG_MAX_POINTS, trianglePruneBins=False)  # make signature factory
        sig_fac.SetBins(THIS_MODULE.FG_BINS)
        sig_fac.Init()
        THIS_MODULE.SIG_FAC = sig_fac

        THIS_MODULE._initialized = True
    elif not THIS_MODULE._initialized:
        delattr(THIS_MODULE, 'INITIALIZED')
        init()

init()
