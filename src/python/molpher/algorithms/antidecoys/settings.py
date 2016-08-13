import os

from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory

from molpher.algorithms.settings import Settings


class AntidecoysSettings(Settings):

    def __init__(
            self
            , source
            , target
            , storage_dir = os.path.abspath('antidecoys_data') # path to a directory where the results will be stored
            , max_threads = None # maximum number of threads to use in parallel computations
            , tree_params = None # custom parameters (same for both trees)
            , max_iters = 100  # maximum number of iterations to spend looking for a single path
            , verbose = False # require verbose output
            , fg_bins = ((0, 2), (2, 5), (5, 8)) # distance bins in the pharmacophore fingerprint
            , fg_min_points = 2 # min number of features encoded in the pharmacophore fingerprint
            , fg_max_points = 3 # max number of features encoded in the pharmacophore fingerprint
            , min_accepted = 1000 # minimum number of morphs the filter will accept on every iteration
            , common_bits_max_thrs = 0.75 # maximum common bits percentage the filter will accept on every iteration
            , common_bits_mean_thrs = 0.5 # if for the mols selected by the filter the mean common bits percentage falls below this value, antidecoys will be turned off
            , antidecoys_min_iters = 10 # minimum number of iterations where antidecoys are optimized
            , antidecoys_max_iters = 50 # maximum number of iterations where antidecoys are optimized
            , distance_thrs = 0.2 # turn antidecoys filter off when the distance between two closest molecules from each tree gets below this value
    ):
        super(AntidecoysSettings, self).__init__(
            source
            , target
            , storage_dir
            , max_threads
            , max_iters
            , tree_params
            , verbose
        )
        self.fg_bins = fg_bins
        self.fg_min_points = fg_min_points
        self.fg_max_points = fg_max_points
        self.min_accepted = min_accepted
        self.common_bits_max_thrs = common_bits_max_thrs
        self.common_bits_mean_thrs = common_bits_mean_thrs
        self.antidecoys_min_iters = antidecoys_min_iters
        self.antidecoys_max_iters = antidecoys_max_iters
        self.distance_thrs = distance_thrs

        # stuff for the pharmacophore fingerprints
        self._fdef_file = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef') # get basic feature definitions
        self._feature_factory = ChemicalFeatures.BuildFeatureFactory(self._fdef_file) # make feature factory
        self.signature_factory = SigFactory(self._feature_factory, minPointCount=self.fg_min_points, maxPointCount=self.fg_max_points, trianglePruneBins=False)  # make signature factory
        self.signature_factory.SetBins(self.fg_bins) # set the distance bins
        self.signature_factory.Init()
