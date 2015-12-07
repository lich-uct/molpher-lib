import molpher
import warnings

from molpher.swig_wrappers.core import TraverseCallback, TraverseOper, MolpherMol


class Callback(TraverseCallback):

    def __init__(self, callback):
        super(self.__class__, self).__init__()
        self.callback = callback

    def processMorph(self, morph):
        self.callback(morph)

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

        self.callback_class = Callback

    @property
    def params(self):
        params = molpher.core.ExplorationParameters(parameters=self.getParams())
        return params.param_dict

    @params.setter
    def params(self, params):
        if params.__class__ == molpher.core.ExplorationParameters \
                or params.__class__ == molpher.wrappers.ExplorationParameters:
            self.setParams(params)
        else:
            new = molpher.core.ExplorationParameters(parameters=self.getParams())
            new.param_dict = params
            self.setParams(new)

    @property
    def leaves(self):
        return self.fetchLeaves()

    @property
    def candidates(self):
        return self.getCandidateMorphs()

    @property
    def candidates_mask(self):
        return self.getCandidateMorphsMask()

    @candidates_mask.setter
    def candidates_mask(self, mask):
        self.setCandidateMorphsMask(mask)

    @property
    def thread_count(self):
        return self.getThreadCount()

    @thread_count.setter
    def thread_count(self, val):
        self.setThreadCount(val)

    def traverse(self, callback, start_mol = None):
        if start_mol and type(start_mol) == MolpherMol:
            TraverseOper(self, self.callback_class(callback), start_mol)()
        elif start_mol:
            mol = self.fetchMol(start_mol)
            TraverseOper(self, self.callback_class(callback), mol)()
        else:
            TraverseOper(self, self.callback_class(callback))()
