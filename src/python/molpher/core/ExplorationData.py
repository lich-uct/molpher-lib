"""
This module houses the :class:`ExplorationData` class:

"""

import molpher
from molpher.core import selectors
from molpher.core.MolpherMol import MolpherMol
from molpher.swig_wrappers.core import FingerprintShortDesc, SimCoeffShortDesc, ChemOperShortDesc


class ExplorationData(molpher.swig_wrappers.core.ExplorationData):
    """
    :param other: an instance of `molpher.swig_wrappers.core.ExplorationData` to wrap with this class
    :type other: `molpher.swig_wrappers.core.ExplorationData`
    :param \*\*kwargs: the morphing parameters to be set (can be incomplete)
    :type \*\*kwargs: `dict`

    .. note:: If both ``other`` and ``**kwargs`` are specified,
        then everything in ``**kwargs`` will be applied *after*
        the instance in ``other`` is wrapped.

    This a specialized version of the `molpher.swig_wrappers.core.ExplorationData` proxy class.
    It implements some additional functionality for ease of use from Python.

    It contains all the information needed to initialize
    an :class:`~molpher.core.ExplorationTree.ExplorationTree` instance.
    Additionally, any tree can be transformed into an instance of this class
    by calling the :meth:`~molpher.core.ExplorationTree.ExplorationTree.asData` method.

    One advantage of this class over the :class:`~molpher.core.ExplorationTree.ExplorationTree`
    is that it allows direct modifications of
    the exploration tree structure. This is especially useful when we want to create
    an initial tree topology before the exploration itself.

    ..  warning:: Note that current implementations of the modification
            methods is experimental and may result
            in undefined behaviour. Therefore, it is only recommended
            to use it as a means of setting morphing parameters
            and spawning tree instances or spawning new trees
            from existing ones without the need to create a snapshot file.

    Because it inherits from `molpher.swig_wrappers.core.ExplorationData`,
    it provides the same interface as the corresponding C++ class,
    but exposes the morphing parameters as object attributes for ease of use.
    These attributes follow a slightly different name convention than the corresponding getters
    and setters of the parent class.
    Their names are derived from the names of the parameters used in the :term:`XML template` files
    that are more self-explanatory and easier to remember and type.
    The table below gives an overview of all available parameters,
    their default values and short descriptions and the respective getters and setters
    of the base class:

    .. include:: param_table.rst

    .. seealso:: `molpher.swig_wrappers.core.ExplorationData`

    """

    class UnknownParameterException(Exception):
        """
        Indicates that an unknown parameter was supplied.
        """

        pass

    def __init__(self, other=None, **kwargs):
        super(ExplorationData, self).__init__()
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

        if other:
            self.this = other.this
        if kwargs:
            self._update_instance(kwargs)

    @staticmethod
    def _parse_options(options):
        """
        A prepossessing step to the :meth:`_update_instance` method.
        Transforms some of the values into values acceptable by
        the setter methods of the base class

        :param options: morphing parameters supplied by caller
        :type options: `dict`
        :return: a transformed dictionary of options
        :rtype: `dict`
        """
        # TODO: add some error checking

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
        """
        Parses a `dict` of :term:`morphing parameters`
        and updates this instance accordingly.

        :param options: morphing parameters to set
        :type options: `dict`
        """

        options = self._parse_options(options)
        for key in options:
            if key in self._SETTERS_MAP:
                self._SETTERS_MAP[key](options[key])
            else:
                raise self.UnknownParameterException('Unknown option: {0}'.format(key))


    @property
    def param_dict(self):
        """
        Holds a dictionary of current :term:`morphing parameters` values for this instance.
        A new dictionary of parameters can be assigned to change them.

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
        The instance becomes invalid, if there are any bad or nonsensical parameter values,
        values are missing (such as undefined :term:`chemical operators`) or the tree structure
        is for any reason unacceptable.

        :return: `True` for a valid instance, `False` for invalid
        :rtype: `bool`
        """

        return self.isValid()

    @property
    def source(self):
        """
        The :term:`source molecule`. All morphs in an :term:`exploration tree` are derived from this
        molecule during morphing. This is the root of the created tree.

        Can be set using a :py:class:`~molpher.core.MolpherMol.MolpherMol` instance
        or a SMILES string of the new :term:`source molecule`.

        :return: current :term:`source molecule`
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
        The :term:`target molecule`. This is the molecule being searched for during morphing.
        In the original version of the algorithm the goal is to
        maximize similarity (minimize structural distance) of the generated morphs and this molecule.

        Can be set using a :py:class:`~molpher.core.MolpherMol.MolpherMol` instance
        or a SMILES string of the new :term:`target molecule`.

        :return: current :term:`target molecule`
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
        A set of :term:`chemical operators` to use. These define how the input molecule and its descendants
        can be manipulated during morphing.

        Can be set using an iterable of the appropriate :term:`selectors` or their names as `str`.
        Any duplicates are automatically removed

        .. include:: oper_table.rst

        :return: names of the current :term:`chemical operators`
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
                # FIXME: check if the correct selectors were supplied
                chosen_selectors.add(selector)
        self._SETTERS_MAP['operators'](tuple(chosen_selectors))

    @property
    def fingerprint(self):
        """
        Returns an identifier of the currently used :term:`molecular fingerprint`.

        .. include:: fing_table.rst

        :return: :term:`molecular fingerprint` identifier
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
        Returns an identifier of the currently used :term:`similarity measure`.

        .. include:: sim_table.rst

        :return: :term:`similarity measure` identifier
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
        If `FilterMorphsOper.WEIGHT` filter is used on an :term:`exploration tree`,
        this will be the minimum weight
        of the :term:`candidate morphs` accepted during a filtering procedure.

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
        If `FilterMorphsOper.WEIGHT` filter is used on an :term:`exploration tree`,
        this will be the maximum weight
        of the :term:`candidate morphs` accepted during a filtering procedure.

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
        The maximum number of morphs allowed to be connected to the tree upon one call to `extend()`.

        If more than `accept_max` morphs with `True` in the appropriate position of `candidates_mask`
        are present in `candidates` and
        `extend()` is called, only first `accept_max` morphs from `candidates` will
        be connected to the tree and the rest will be discarded.

        .. seealso:: `ExplorationTree.extend()`

        :return: maximum number of candidates accepted upon `extend()`
        :rtype: `int`
        """

        return self._GETTERS_MAP['accept_max']()

    @accept_max.setter
    def accept_max(self, value):
        self._SETTERS_MAP['accept_max'](value)

    @property
    def far_produce(self):
        """
        The maximum number of morphs generated from one leaf when the leaf of the tree currently being
        processed with `generateMorphs()` lies more than `far_close_threshold` from
        the :term:`target molecule`.

        .. seealso:: `ExplorationTree.generateMorphs()`

        :return: maximum number of morphs to produce with a `generateMorphs()` call
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
        processed with `generateMorphs()` lies less than `far_close_threshold` from
        the :term:`target molecule`.

        .. seealso:: `ExplorationTree.generateMorphs()`

        :return: maximum number of morphs to produce with an `generateMorphs()` call
        :rtype: `int`
        """

        return self._GETTERS_MAP['close_produce']()

    @close_produce.setter
    def close_produce(self, value):
        self._SETTERS_MAP['close_produce'](value)

    @property
    def far_close_threshold(self):
        """
        This distance threshold controls the number of :term:`morphs <morph>` generated
        with `generateMorphs()` for molecules closer
        or further from the :term:`target molecule`. :term:`Morphs <morph>` that
        have distance from the :term:`target molecule` lower than `far_close_threshold`
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
        descendants and none of them are closer to the :term:`target molecule` than the molecule in question, then
        the molecule is permanently removed from the tree with all of its descendants when `prune()`
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
        A molecule that has not produced any morphs closer to the :term:`target molecule` than itself
        (a :term:`non-producing molecule`) for `non_producing_survive` number of calls to `extend()`
        will have its descendants removed during the next `prune()` call.

        .. seealso:: `MolpherMol.getItersWithoutDistImprovement()`

        :return: number of calls  to `generateMorphs()` before descendants of
            a :term:`non-producing molecule`
            are removed from the tree
        :rtype: `int`
        """

        return self._GETTERS_MAP['non_producing_survive']()

    @non_producing_survive.setter
    def non_producing_survive(self, value):
        self._SETTERS_MAP['non_producing_survive'](value)

    @staticmethod
    def load(snapshot):
        """
        A factory method to create an instance of :class:`ExplorationData`
        from a :term:`tree snapshot`.

        :param snapshot: path to the snapshot file
        :type snapshot: `str`
        :return: new instance representing the data loaded from the snapshot file
        :rtype: :class:`ExplorationData`
        """

        data = super(ExplorationData, ExplorationData).load(snapshot)
        return ExplorationData(other=data)