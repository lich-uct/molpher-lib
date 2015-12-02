import molpher
import warnings

class ExplorationTree(molpher.swig_wrappers.core.ExplorationTree):

    def __init__(self, params=None, source=None, target=None):
        if params and (source or target):
            warnings.warn(
                "Both parameters and a source or a target specified. Using the values in parameters..."
                , RuntimeWarning
            )
        if params and (params.__class__ == molpher.core.ExplorationParameters
                or params.__class__ == molpher.wrappers.ExplorationParameters):
            super(ExplorationTree, self).__init__(params)
        elif params:
            _params = molpher.core.ExplorationParameters(**params)
            super(ExplorationTree, self).__init__(_params)
        elif source:
            if not target:
                super(ExplorationTree, self).__init__(source)
            else:
                raise NotImplementedError('ExplorationTree constructor does not '
                                          'support explicit target initialization, yet. '
                                          'Use an ExplorationParameters instance instead.'
                                          )
        else:
            raise RuntimeError('You must specify either `params` or `source`.')

    @property
    def params(self):
        params = self.getParams()
        return {
            'source' : params.getSourceMol().getSMILES()
            , 'target' : params.getTargetMol().getSMILES()
            , 'operators' : params.getChemOperators()
            , 'fingerprint' : params.getFingerprint()
            , 'similarity' : params.getSimilarityCoef()
            , 'weight_min' : params.getMinAcceptableMolecularWeight()
            , 'weight_max' : params.getMaxAcceptableMolecularWeight()
            , 'accept_min' : params.getCntCandidatesToKeep()
            , 'accept_max' : params.getCntCandidatesToKeepMax()
            , 'far_produce' : params.getCntMorphs()
            , 'close_produce' : params.getCntMorphsInDepth()
            , 'far_close_threshold' : params.getDistToTargetDepthSwitch()
            , 'max_morphs_total' : params.getCntMaxMorphs()
            , 'non_producing_survive' : params.getItThreshold()
        }

    @params.setter
    def params(self, params):
        if params.__class__ == molpher.core.ExplorationParameters \
                or params.__class__ == molpher.wrappers.ExplorationParameters:
            self.setParams(params)
        else:
            _params = molpher.core.ExplorationParameters(**params)
            self.setParams(_params)