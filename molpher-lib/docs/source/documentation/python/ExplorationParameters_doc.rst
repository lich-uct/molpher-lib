ExplorationParameters
=====================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExplorationParameters

    Contains information about the `exploration parameters` and provides facilities to read them
    and change them. The constructor does not take any arguments and automatically initilizes
    the parameters with the default values (see `molpher.core.ExplorationParameters` for details).

    .. automethod:: getSourceMol

        Returns the `source molecule`.

        :return: `MolpherMol` instance representing the current `source molecule`
        :rtype: `MolpherMol`

    .. automethod:: getTargetMol

        Returns the `target molecule`.

        :return: `MolpherMol` instance representing the current `target molecule`
        :rtype: `MolpherMol`

    .. automethod:: getChemOperators

        Returns a `tuple` of current `chemical operators` as `str` identifiers.

        .. todo:: finish the table

        .. include:: oper_table.rst

        :return: current chemical operators
        :rtype: `tuple` of `str`

    .. automethod:: getFingerprint

        Returns an identifier of the currently set `molecular fingerprint`.

        .. todo:: add table

        :return: `molecular fingerprint` identifier
        :rtype: `str`

    .. automethod:: getSimilarityCoef

        Returns an identifier of the currently set `similarity measure`.

        .. todo:: add table

        :return: `similarity measure` identifier
        :rtype: `str`

    .. automethod:: getMinAcceptableMolecularWeight

        If `FilterMorphsOper.WEIGHT` filter is used on an `exploration tree` with the current parameters set,
        this will be the minimum weight
        of the `candidate morphs` accepted during a filtering procedure.

        .. seealso:: `ExplorationTree.filterMorphs()`

        .. note:: equivalent to ``weight_min`` in `molpher.core.ExplorationParameters`

        :return: minimum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getMaxAcceptableMolecularWeight

        If `FilterMorphsOper.WEIGHT` filter is used on an `exploration tree` with the current parameters set,
        this will be the maximum weight
        of the `candidate morphs` accepted during a filtering procedure.

        .. seealso:: `ExplorationTree.filterMorphs()`

        .. note:: equivalent to ``weight_max`` in `molpher.core.ExplorationParameters`

        :return: maximum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getCntCandidatesToKeep

        .. todo:: document

        .. note:: equivalent to ``accept_min`` in `molpher.core.ExplorationParameters`

        :return: min number of candidates
        :rtype: `int`

    .. automethod:: getCntCandidatesToKeepMax

        .. todo:: document

        .. note:: equivalent to ``accept_max`` in `molpher.core.ExplorationParameters`

        :return: max number of candidates
        :rtype: `int`

    .. automethod:: getCntMorphs

        .. todo:: document

        .. note:: equivalent to ``far_produce`` in `molpher.core.ExplorationParameters`

        :return:
        :rtype: `int`

    .. automethod:: getCntMorphsInDepth

        .. todo:: document

        .. note:: equivalent to ``close_produce`` in `molpher.core.ExplorationParameters`

        :return:
        :rtype: `int`

    .. automethod:: getDistToTargetDepthSwitch

        .. todo:: document

        .. note:: equivalent to ``far_close_threshold`` in `molpher.core.ExplorationParameters`

        :return:
        :rtype: `float`

    .. automethod:: getCntMaxMorphs

        .. todo:: document

        .. note:: equivalent to ``max_morhps_total`` in `molpher.core.ExplorationParameters`

        :return:
        :rtype: `int`

    .. automethod:: getItThreshold

        .. todo:: document

        .. note:: equivalent to ``non_producing_survive`` in `molpher.core.ExplorationParameters`

        :return:
        :rtype: `int`