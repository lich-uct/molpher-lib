ExplorationData
=====================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExplorationData

    :param \*args: either another instance of this class or empty

    .. note::  An instance of this class can be either initilized automatically (no arguments
        provided to the constructor) or from another instance (passed as a single constructor parameter).
        If the class is initilized automatically (no parameters passed to the constructor),
        the `exploration parameters` are set to the default values (see the description of `molpher.core.ExplorationData`
        for an overview of default parameter values).

    Contains information about the `exploration parameters` and provides facilities to read them
    and change them.

    This class uses getter and setter methods to read and modify the parameter values. The
    naming conventions tightly follow the wrapped C++ implementation and should be consistent
    with the parameter names published in [1]_.

    .. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

    .. automethod:: getSource

        Returns the `source molecule`.

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.source`

        :return: `MolpherMol` instance representing the current `source molecule`
        :rtype: `MolpherMol`

    .. automethod:: getTarget

        Returns the `target molecule`.

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.target`

        :return: `MolpherMol` instance representing the current `target molecule`
        :rtype: `MolpherMol`

    .. automethod:: getChemicalOperators

        Returns a `tuple` of current `chemical operators` as `str` identifiers.

        .. include:: oper_table.rst

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.operators`

        :return: current chemical operators
        :rtype: `tuple` of `str`

    .. automethod:: getFingerprint

        Returns an identifier of the currently set `molecular fingerprint`.

        .. include:: fing_table.rst

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.fingerprint`

        :return: `molecular fingerprint` identifier
        :rtype: `str`

    .. automethod:: getSimilarityCoefficient

        Returns an identifier of the currently set `similarity measure`.

        .. include:: sim_table.rst

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.similarity`

        :return: `similarity measure` identifier
        :rtype: `str`

    .. automethod:: getMinAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.weight_min`

        :return: minimum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getMaxAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.weight_max`

        :return: maximum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getCntCandidatesToKeep

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.accept_min`

        :return: min number of candidates
        :rtype: `int`

    .. automethod:: getCntCandidatesToKeepMax

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.accept_max`

        :return: max number of candidates
        :rtype: `int`

    .. automethod:: getCntMorphs

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.far_produce`

        :return:
        :rtype: `int`

    .. automethod:: getCntMorphsInDepth

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.close_produce`

        :return: number of morphs to generate
        :rtype: `int`

    .. automethod:: getDistToTargetDepthSwitch

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.far_close_threshold`

        :return: distance threshold
        :rtype: `float`

    .. automethod:: getCntMaxMorphs

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.max_morphs_total`

        :return: maximum number of non-producing morphs per molecule
        :rtype: `int`

    .. automethod:: getItThreshold

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.non_producing_survive`

        :return: number of generations until non-producing morphs are pruned
        :rtype: `int`

    .. automethod:: setSource

        Set the `source molecule`.

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.source`

        :param \*args: `MolpherMol` instance or SMILES representing the current `source molecule`
        :type \*args: `MolpherMol` or `str`

    .. automethod:: setTarget

        Set the `target molecule`.

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.target`

        :param \*args: `MolpherMol` instance or SMILES representing the current `source molecule`
        :type \*args: `MolpherMol` or `str`

    .. automethod:: setChemicalOperators

        Set the current `chemical operators` as `str` identifiers.

        .. include:: oper_table.rst

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.operators`

        :param \*args: current `chemical operators` as a `tuple` of `str` identifiers
        :type \*args: `tuple` of `str`

    .. automethod:: setFingerprint

        Set the `molecular fingerprint`.

        .. include:: fing_table.rst

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.fingerprint`

        :param \*args: `molecular fingerprint` identifier
        :type \*args: `str`

    .. automethod:: setSimilarityCoefficient

        Set the `similarity measure`.

        .. include:: sim_table.rst

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.similarity`

        :param \*args: `similarity measure` identifier
        :type \*args: `str`

    .. automethod:: setMinAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.weight_min`

        :param \*args: minimum acceptable weight during filtering
        :type \*args: `float`

    .. automethod:: setMaxAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.weight_max`

        :param \*args: maximum acceptable weight during filtering
        :type \*args: `float`

    .. automethod:: setCntCandidatesToKeep

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.accept_min`

        :param \*args: min number of candidates
        :type \*args: `int`

    .. automethod:: setCntCandidatesToKeepMax

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.accept_max`

        :param \*args: max number of candidates
        :type \*args: `int`

    .. automethod:: setCntMorphs

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.far_produce`

        :param \*args:
        :type \*args: `int`

    .. automethod:: setCntMorphsInDepth

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.close_produce`

        :param \*args: number of morphs to generate
        :type \*args: `int`

    .. automethod:: setDistToTargetDepthSwitch

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.far_close_threshold`

        :param \*args: distance threshold
        :type \*args: `float`

    .. automethod:: setCntMaxMorphs

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.max_morphs_total`

        :param \*args: maximum number of non-producing morphs per molecule
        :type \*args: `int`

    .. automethod:: setItThreshold

        .. seealso:: `molpher.core.ExplorationData.ExplorationData.non_producing_survive`

        :param \*args: number of generations until non-producing morphs are pruned
        :type \*args: `int`