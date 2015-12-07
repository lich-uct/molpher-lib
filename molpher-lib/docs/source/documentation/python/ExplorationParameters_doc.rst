ExplorationParameters
=====================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExplorationParameters

    :param \*args: either another instance of this class or empty

    Contains information about the `exploration parameters` and provides facilities to read them
    and change them.

    An instance of this class can be either initilized automatically (no arguments
    provided to the constructor) or from another instance (passed as a single constructor parameter).
    If the class is initilized automatically (no parameters passed to the constructor),
    the `exploration parameters` are set to the default values (see the description of `molpher.core.ExplorationParameters`
    for an overview of default parameters).

    .. automethod:: getSourceMol

        Returns the `source molecule`.

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.source`

        :return: `MolpherMol` instance representing the current `source molecule`
        :rtype: `MolpherMol`

    .. automethod:: getTargetMol

        Returns the `target molecule`.

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.target`

        :return: `MolpherMol` instance representing the current `target molecule`
        :rtype: `MolpherMol`

    .. automethod:: getChemOperators

        Returns a `tuple` of current `chemical operators` as `str` identifiers.

        .. todo:: finish the table

        .. include:: oper_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.operators`

        :return: current chemical operators
        :rtype: `tuple` of `str`

    .. automethod:: getFingerprint

        Returns an identifier of the currently set `molecular fingerprint`.

        .. include:: fing_table.rst

        .. todo:: finish table

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.fingerprint`

        :return: `molecular fingerprint` identifier
        :rtype: `str`

    .. automethod:: getSimilarityCoef

        Returns an identifier of the currently set `similarity measure`.

        .. include:: sim_table.rst

        .. todo:: finish table

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.similarity`

        :return: `similarity measure` identifier
        :rtype: `str`

    .. automethod:: getMinAcceptableMolecularWeight

        If `FilterMorphsOper.WEIGHT` filter is used on an `exploration tree` with the current parameters set,
        this will be the minimum weight
        of the `candidate morphs` accepted during a filtering procedure.

        .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.weight_min`

        :return: minimum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getMaxAcceptableMolecularWeight

        If `FilterMorphsOper.WEIGHT` filter is used on an `exploration tree` with the current parameters set,
        this will be the maximum weight
        of the `candidate morphs` accepted during a filtering procedure.

        .. seealso:: .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.weight_max`

        :return: maximum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getCntCandidatesToKeep

        .. todo:: document

        .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.accept_min`

        :return: min number of candidates
        :rtype: `int`

    .. automethod:: getCntCandidatesToKeepMax

        .. todo:: document

        .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.accept_max`

        :return: max number of candidates
        :rtype: `int`

    .. automethod:: getCntMorphs

        .. todo:: document

        .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.far_produce`

        :return:
        :rtype: `int`

    .. automethod:: getCntMorphsInDepth

        .. todo:: document

        .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.close_produce`

        :return:
        :rtype: `int`

    .. automethod:: getDistToTargetDepthSwitch

        .. todo:: document

        .. seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationParameters.ExplorationParameters.far_close_threshold`

        :return:
        :rtype: `float`

    .. automethod:: getCntMaxMorphs

        .. todo:: document

        .. seealso:: `ExplorationTree.prune()` and `molpher.core.ExplorationParameters.ExplorationParameters.max_morphs_total`

        :return:
        :rtype: `int`

    .. automethod:: getItThreshold

        .. todo:: document

        .. seealso:: `ExplorationTree.prune()` and `molpher.core.ExplorationParameters.ExplorationParameters.non_producing_survive`

        :return:
        :rtype: `int`