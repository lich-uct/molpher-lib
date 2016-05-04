import molpher
from molpher.core import selectors
from molpher.core.MolpherMol import MolpherMol
from molpher.swig_wrappers.core import FingerprintShortDesc, SimCoeffShortDesc, ChemOperShortDesc


class ExplorationData(molpher.swig_wrappers.core.ExplorationData):
    """
    :param parameters: another instance to copy the parameters from
    :type parameters: `molpher.swig_wrappers.core.ExplorationData` or its derived class
    :param \*\*kwargs: morphing parameters
    :type \*\*kwargs: `dict`

    .. note:: If both ``parameters`` and ``**kwargs`` are specified while intilizing an instance of this class,
        the parameters are initilized from the instance
        provided in ``parameters``, but with the parameters in ``**kwargs`` taking precedence.

    This class houses morphing parameters and can be used to initialize a `molpher.core.ExplorationTree` instance.

    It inherits from `molpher.swig_wrappers.core.ExplorationData`. Therefore, it provides the same interface.
    However, it also exposes the parameters as object attributes that use a slightly different name convention
    (derived from the names of the parameters used in the `XML template` files) that should make the
    parameter names more self-explanatory.

    The table below gives an overview of all available parameters, their respective getters and setters
    in the base class, default values and short descriptions. If you want to know more, refer to
    the individual attribute descriptions below.

    .. include:: param_table.rst

    """

    class UnknownParameterException(Exception):
        """
        Indicates that an unknown parameter was supplied.
        """

        pass

    def __init__(self, **kwargs):
        super(ExplorationData, self).__init__()
        other = None
        if 'other' in kwargs:
            other = kwargs.pop('other')
        self._SETTERS_MAP = {
            'source' : self.setSource
            , 'target' : self.setTarget
            , 'operators' : self.setChemicalOperators
            , 'fingerprint' : self.setFingerprint
            , 'similarity' : self.setSimilarityCoefficient
            , 'weight_min' : self.setMinAcceptableMolecularWeight
            , 'weight_max' : self.setMaxAcceptableMolecularWeight
            , 'accept_min' : self.setCntCandidatesToKeep
            , 'accept_max' : self.setCntCandidatesToKeepMax
            , 'far_produce' : self.setCntMorphs
            , 'close_produce' : self.setCntMorphsInDepth
            , 'far_close_threshold' : self.setDistToTargetDepthSwitch
            , 'max_morphs_total' : self.setCntMaxMorphs
            , 'non_producing_survive' : self.setItThreshold
        }
        self._GETTERS_MAP = {
            'source' : self.getSource
            , 'target' : self.getTarget
            , 'operators' : self.getChemicalOperators
            , 'fingerprint' : self.getFingerprint
            , 'similarity' : self.getSimilarityCoefficient
            , 'weight_min' : self.getMinAcceptableMolecularWeight
            , 'weight_max' : self.getMaxAcceptableMolecularWeight
            , 'accept_min' : self.getCntCandidatesToKeep
            , 'accept_max' : self.getCntCandidatesToKeepMax
            , 'far_produce' : self.getCntMorphs
            , 'close_produce' : self.getCntMorphsInDepth
            , 'far_close_threshold' : self.getDistToTargetDepthSwitch
            , 'max_morphs_total' : self.getCntMaxMorphs
            , 'non_producing_survive' : self.getItThreshold
        }

        if not other:
            self._update_instance(kwargs)
        else:
            self.this = other.this

    @staticmethod
    def _parse_options(options):
        check = lambda key : key in options and type(options[key]) == str
        check_iterable = lambda key : key in options and [x for x in options[key] if type(x) == str]
        if options:
            if check('source'):
                options['source'] = MolpherMol(options['source'])
            if check('target'):
                options['target'] = MolpherMol(options['target'])
            if check_iterable('operators'):
                options['operators'] = tuple( getattr(selectors, operator) for operator in options['operators'] )
            if check('fingerprint'):
                options['fingerprint'] = getattr(selectors, options['fingerprint'])
            if check('similarity'):
                options['similarity'] = getattr(selectors, options['similarity'])
        return options


    def _update_instance(self, options):
        options = self._parse_options(options)
        for key in options:
            if key in self._SETTERS_MAP:
                self._SETTERS_MAP[key](options[key])
            else:
                raise self.UnknownParameterException('Unknown option: {0}'.format(key))


    @property
    def param_dict(self):
        """
        Holds a dictionary of current parameter values for this instance.
        A new dictionary of parameters can be assigned to change the current parameters.

        :return: a dictionary of parameters
        :rtype: `dict`
        """

        return {
            'source' : self._GETTERS_MAP['source']().getSMILES()
            , 'target' : self._GETTERS_MAP['target']().getSMILES()
            , 'operators' : tuple( ChemOperShortDesc(operator) for operator in self._GETTERS_MAP['operators']() )
            , 'fingerprint' : FingerprintShortDesc(self._GETTERS_MAP['fingerprint']())
            , 'similarity' : SimCoeffShortDesc(self._GETTERS_MAP['similarity']())
            , 'weight_min' : self._GETTERS_MAP['weight_min']()
            , 'weight_max' : self._GETTERS_MAP['weight_max']()
            , 'accept_min' : self._GETTERS_MAP['accept_min']()
            , 'accept_max' : self._GETTERS_MAP['accept_max']()
            , 'far_produce' : self._GETTERS_MAP['far_produce']()
            , 'close_produce' : self._GETTERS_MAP['close_produce']()
            , 'far_close_threshold' : self._GETTERS_MAP['far_close_threshold']()
            , 'max_morphs_total' : self._GETTERS_MAP['max_morphs_total']()
            , 'non_producing_survive' : self._GETTERS_MAP['non_producing_survive']()
        }

    @param_dict.setter
    def param_dict(self, options):
        self._update_instance(options)

    @property
    def is_valid(self):
        """
        Shows if this instance represents valid parameters.
        The instance becomes invalid, if there are any bad or nonsensical parameter values
        or some values are missing (such as undefined `chemical operators`).

        :return: whether validity is `True` or `False`
        :rtype: `bool`
        """

        return self.isValid()

    @property
    def source(self):
        """
        The `source molecule`. All morphs in an `exploration tree` are derived from this
        molecule during morphing. This is the root of the tree.

        :return: :py:class:`~molpher.core.MolpherMol.MolpherMol` instance or SMILES string representing the current `source molecule`
        :rtype: :py:class:`~molpher.core.MolpherMol.MolpherMol`
        """

        return MolpherMol(other=self._GETTERS_MAP['source']())

    @source.setter
    def source(self, value):
        if type(value) == str:
            mol = MolpherMol(value)
            self._SETTERS_MAP['source'](mol)
        elif isinstance(value, MolpherMol):
            self._SETTERS_MAP['source'](value)
        else:
            raise Exception("Invalid input. Need 'str' or 'MolpherMol'...")

    @property
    def target(self):
        """
        The `target molecule`. This is the molecule being searched for during morphing. The goal is to
        maximize similarity (minimize distance) of the generated morphs and this molecule.

        :return: :py:class:`~molpher.core.MolpherMol.MolpherMol` instance representing the current `target molecule`
        :rtype: :py:class:`~molpher.core.MolpherMol.MolpherMol`
        """

        return MolpherMol(other=self._GETTERS_MAP['target']())

    @target.setter
    def target(self, value):
        if type(value) == str:
            mol = MolpherMol(value)
            self._SETTERS_MAP['target'](mol)
        elif isinstance(value, MolpherMol):
            self._SETTERS_MAP['target'](value)
        else:
            raise Exception("Invalid input. Need 'str' or 'MolpherMol'...")

    @property
    def operators(self):
        """

        A set of `chemical operators` to use. These define how the input molecule and its descendents
        can be manipulated during morphing.

        .. include:: oper_table.rst

        :return: current chemical operators
        :rtype: `tuple` of `str`
        """

        return tuple( ChemOperShortDesc(operator) for operator in self._GETTERS_MAP['operators']() )

    @operators.setter
    def operators(self, value):
        chosen_selectors = set()
        for selector in value:
            if type(selector) == str:
                chosen_selectors.add(getattr(selectors, selector))
            else:
                chosen_selectors.add(selector)
        self._SETTERS_MAP['operators'](tuple(chosen_selectors))

    @property
    def fingerprint(self):
        """
        Returns an identifier of the currently set `molecular fingerprint`.

        .. include:: fing_table.rst

        :return: `molecular fingerprint` identifier
        :rtype: `str`
        """

        return FingerprintShortDesc(self._GETTERS_MAP['fingerprint']())

    @fingerprint.setter
    def fingerprint(self, value):
        if type(value) == str:
            value = getattr(selectors, value)
        self._SETTERS_MAP['fingerprint'](value)

    @property
    def similarity(self):
        """
        Returns an identifier of the currently set `similarity measure`.

        .. include:: sim_table.rst

        :return: `similarity measure` identifier
        :rtype: `str`
        """

        return SimCoeffShortDesc(self._GETTERS_MAP['similarity']())

    @similarity.setter
    def similarity(self, value):
        if type(value) == str:
            value = getattr(selectors, value)
        self._SETTERS_MAP['similarity'](value)

    @property
    def weight_min(self):
        """
        If `FilterMorphsOper.WEIGHT` filter is used on an `exploration tree` with the current parameters set,
        this will be the minimum weight
        of the `candidate morphs` accepted during a filtering procedure.

        .. seealso:: `ExplorationTree.filterMorphs()`

        :return: minimum acceptable weight during filtering
        :rtype: `float`
        """

        return self._GETTERS_MAP['weight_min']()

    @weight_min.setter
    def weight_min(self, value):
        self._SETTERS_MAP['weight_min'](value)

    @property
    def weight_max(self):
        """
        If `FilterMorphsOper.WEIGHT` filter is used on an `exploration tree` with the current parameters set,
        this will be the maximum weight
        of the `candidate morphs` accepted during a filtering procedure.

        .. seealso:: `ExplorationTree.filterMorphs()`

        :return: maximum acceptable weight during filtering
        :rtype: `float`
        """

        return self._GETTERS_MAP['weight_max']()

    @weight_max.setter
    def weight_max(self, value):
        self._SETTERS_MAP['weight_max'](value)

    @property
    def accept_min(self):
        """
        If `FilterMorphsOper.PROBABILITY` is used during filtering, this is the number of morphs
        accepted with 100% probability.

        .. seealso:: `ExplorationTree.filterMorphs()`

        :return: minimum number of candidates accepted during probability filtering
        :rtype: `int`
        """

        return self._GETTERS_MAP['accept_min']()

    @accept_min.setter
    def accept_min(self, value):
        self._SETTERS_MAP['accept_min'](value)

    @property
    def accept_max(self):
        """
        This is the maximum number of morphs allowed to be connected to the tree.
        If more than `accept_max` morphs with `True` on the appropriate position in `ExplorationTree.candidates_mask`
        are present in `ExplorationTree.candidates` and
        `ExplorationTree.extend()` is called, only first `accept_max` morphs from `ExplorationTree.candidates` will
        be connected to the tree and the rest will be discarded.

        .. seealso:: `ExplorationTree.extend()`

        :return: maximum number of candidates accepted upon `ExplorationTree.extend()`
        :rtype: `int`
        """

        return self._GETTERS_MAP['accept_max']()

    @accept_max.setter
    def accept_max(self, value):
        self._SETTERS_MAP['accept_max'](value)

    @property
    def far_produce(self):
        """
        This is the maximum number of morphs generated from one leaf when the leaf of the tree currently being
        processed with `ExplorationTree.generateMorphs()` lies more than `far_close_threshold` from
        the `target molecule`.

        :return: maximum number of morphs to produce with an `ExplorationTree.generateMorphs()` call
        :rtype: `int`
        """

        return self._GETTERS_MAP['far_produce']()

    @far_produce.setter
    def far_produce(self, value):
        self._SETTERS_MAP['far_produce'](value)

    @property
    def close_produce(self):
        """
        This is the maximum number of morphs generated from one leaf when the leaf of the tree currently being
        processed with `ExplorationTree.generateMorphs()` lies less than `far_close_threshold` from
        the `target molecule`.

        :return: maximum number of morphs to produce with an `ExplorationTree.generateMorphs()` call
        :rtype: `int`
        """

        return self._GETTERS_MAP['close_produce']()

    @close_produce.setter
    def close_produce(self, value):
        self._SETTERS_MAP['close_produce'](value)

    @property
    def far_close_threshold(self):
        """
        This distance threshold controls the number of `morphs <morph>` generated
        with `ExplorationTree.generateMorphs()` for molecules closer
        or further from the `target molecule`. `Morphs <morph>` that
        have distance from the `target molecule` lower than `far_close_threshold`
        are considered to be close.

        .. seealso:: `far_produce` and `close_produce`

        :return: distance threshold for `far_produce` and `close_produce`
        :rtype: `float`
        """

        return self._GETTERS_MAP['far_close_threshold']()

    @far_close_threshold.setter
    def far_close_threshold(self, value):
        self._SETTERS_MAP['far_close_threshold'](value)

    @property
    def max_morphs_total(self):
        """
        This value is the maximum number of morphs allowed to be generated from one molecule.
        If the number of generated morphs exceeds this number, all additional morphs can be filtered
        out using the `FilterMorphsOper.MAX_DERIVATIONS` filter.

        It is also the maximum number of 'bad morphs' generated from one molecule. If a molecule has more than `max_morphs_total`
        descendents and none of them are closer to the `target molecule` than the molecule in question, then
        the molecule is permanently removed from the tree with all of its descendents when `ExplorationTree.prune()`
        is called.

        .. seealso:: `ExplorationTree.filterMorphs()` and `ExplorationTree.prune()`

        :return: maximum number of 'bad morphs' before pruning
        :rtype: `int`
        """

        return self._GETTERS_MAP['max_morphs_total']()

    @max_morphs_total.setter
    def max_morphs_total(self, value):
        self._SETTERS_MAP['max_morphs_total'](value)

    @property
    def non_producing_survive(self):
        """
        A molecule that has not produced any morphs closer to the `target molecule` than itself (a `non-producing molecule`)
        for `non_producing_survive` number of calls to `ExplorationTree.extend()`
        will have its descendents removed during the next `ExplorationTree.prune()` call.

        .. seealso:: `MolpherMol.getItersWithoutDistImprovement()`

        :return: number of calls  to `ExplorationTree.generateMorphs()` before descendents of a
            'non-producing molecule'
            are removed from the tree
        :rtype: `int`
        """

        return self._GETTERS_MAP['non_producing_survive']()

    @non_producing_survive.setter
    def non_producing_survive(self, value):
        self._SETTERS_MAP['non_producing_survive'](value)

    @staticmethod
    def load(snapshot):
        data = super(ExplorationData, ExplorationData).load(snapshot)
        return ExplorationData(other=data)