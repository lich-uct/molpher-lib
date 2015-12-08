import molpher
import warnings

from molpher.swig_wrappers.core import TraverseCallback, TraverseOper, MolpherMol


class Callback(TraverseCallback):
    """
    :param callback: callable with one positional argument
    :type callback: any callable object

    Basic callback class used to traverse the tree in `ExplorationTree.traverse()` method.

    It registers a callable and calls it every time a `morph` is processed.

    """

    def __init__(self, callback):
        super(self.__class__, self).__init__()
        self._callback = callback

    def processMorph(self, morph):
        """
        Calls the registered callback for a `morph`

        :param morph: morph in a tree
        :type morph: `MolpherMol`
        """

        self._callback(morph)

class ExplorationTree(molpher.swig_wrappers.core.ExplorationTree):
    """
    :param params: the morphing parameters (optional if ``source`` is specified)
    :type params: `molpher.swig_wrappers.core.ExplorationParameters` or any of its derived classes
    :param source: SMILES of the source molecule
    :type source: `str`
    :param target: SMILES of the target molecule (optional)
    :type target: `str`

    ..  note:: When ``params`` are specified, ``source`` and ``target`` are ignored.

    ..  note:: If the target molecule is missing a single carbon atom is supplied instead.

    This specialized version of `molpher.swig_wrappers.core.ExplorationTree`
    implements some additional functionality for ease of use from Python.

    """

    def __init__(self, params=None, source=None, target=None):
        if params and (source or target):
            warnings.warn(
                "Both parameters and a source or a target specified. Using the values in parameters..."
                , RuntimeWarning
            )
        if params and (params.__class__ == molpher.core.ExplorationParameters
                or params.__class__ == molpher.wrappers.ExplorationParameters): # FIXME: check only for the base class molpher.wrappers.ExplorationParameters
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

        self._callback_class = Callback

    @property
    def params(self):
        """
        A dictionary representing the current `exploration parameters`.

        It is possible to assign a new dictionary to update the current parameters.
        Only parameters defined in the supplied dictionary are changed.

        :return: current parameters
        :rtype: `dict`
        """
        params = molpher.core.ExplorationParameters(parameters=self.getParams())
        return params.param_dict

    @params.setter
    def params(self, params):
        if params.__class__ == molpher.core.ExplorationParameters \
                or params.__class__ == molpher.wrappers.ExplorationParameters: # FIXME: check only for the base class molpher.wrappers.ExplorationParameters
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
            TraverseOper(self, self._callback_class(callback), start_mol)()
        elif start_mol:
            mol = self.fetchMol(start_mol)
            TraverseOper(self, self._callback_class(callback), mol)()
        else:
            TraverseOper(self, self._callback_class(callback))()
