import os

from molpher import random

from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory


class AntidecoysSettings:

    def __init__(
            self
            , seed = None # random seed for molpher, if it needs to be constant
            , storage_dir = os.path.abspath('antidecoys_data') # path to the directory with the results
            , max_threads = 4 # maximum number of threads to use in parallel computations
            , fg_bins = ((0, 2), (2, 5), (5, 8)) # distance bins in the pharmacophore fingerprint
            , fg_min_points = 2 # min number of features encoded in the pharmacophore fingerprint
            , fg_max_points = 3 # max number of features encoded in the pharmacophore fingerprint
            , min_accepted = 1000 # minimum number of morphs the filter will accept on every iteration
            , common_bits_max_thrs = 0.75 # maximum common bits percentage the filter will accept on every iteration
            , common_bits_mean_thrs = 0.5 # if for the mols selected by the filter the mean common bits percentage falls below this value, antidecoys will be turned off
            , min_iters = 10 # minimum number of iterations where antidecoys are optimized
            , max_iters = 50 # maximum number of iterations where antidecoys are optimized
            , distance_thrs = 0.2 # turn antidecoys filter off when the distance between two closest molecules from each tree gets below this value
            , max_iters_per_path = 100  # maximum number of iterations to spend looking for a path
            , max_paths = 50
    ):
        self.seed = seed
        if self.seed not in (False, None):
            random.set_random_seed(self.seed)
        self.storage_dir = storage_dir
        self.max_threads = max_threads
        self.fg_bins = fg_bins
        self.fg_min_points = fg_min_points
        self.fg_max_points = fg_max_points
        self.min_accepted = min_accepted
        self.common_bits_max_thrs = common_bits_max_thrs
        self.common_bits_mean_thrs = common_bits_mean_thrs
        self.min_iters = min_iters
        self.max_iters = max_iters
        self.distance_thrs = distance_thrs
        self.max_iters_per_path = max_iters_per_path
        self.max_paths = max_paths

        # stuff for the pharmacophore fingerprints
        self._fdef_file = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef') # get basic feature definitions
        self._feature_factory = ChemicalFeatures.BuildFeatureFactory(self._fdef_file) # make feature factory
        self.signature_factory = SigFactory(self._feature_factory, minPointCount=self.fg_min_points, maxPointCount=self.fg_max_points, trianglePruneBins=False)  # make signature factory
        self.signature_factory.SetBins(self.fg_bins) # set the distance bins
        self.signature_factory.Init()
