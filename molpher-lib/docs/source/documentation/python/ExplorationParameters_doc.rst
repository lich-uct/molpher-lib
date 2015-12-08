ExplorationParameters
=====================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExplorationParameters

    :param \*args: either another instance of this class or empty

    .. note::  An instance of this class can be either initilized automatically (no arguments
        provided to the constructor) or from another instance (passed as a single constructor parameter).
        If the class is initilized automatically (no parameters passed to the constructor),
        the `exploration parameters` are set to the default values (see the description of `molpher.core.ExplorationParameters`
        for an overview of default parameter values).

    Contains information about the `exploration parameters` and provides facilities to read them
    and change them.

    This class uses getter and setter methods to read and modify the parameter values. The
    naming conventions tightly follow the wrapped C++ implementation and should be consistent
    with the parameter names published in [1]_.

    .. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

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

        .. include:: oper_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.operators`

        :return: current chemical operators
        :rtype: `tuple` of `str`

    .. automethod:: getFingerprint

        Returns an identifier of the currently set `molecular fingerprint`.

        .. include:: fing_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.fingerprint`

        :return: `molecular fingerprint` identifier
        :rtype: `str`

    .. automethod:: getSimilarityCoef

        Returns an identifier of the currently set `similarity measure`.

        .. include:: sim_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.similarity`

        :return: `similarity measure` identifier
        :rtype: `str`

    .. automethod:: getMinAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.weight_min`

        :return: minimum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getMaxAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.weight_max`

        :return: maximum acceptable weight during filtering
        :rtype: `float`

    .. automethod:: getCntCandidatesToKeep

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.accept_min`

        :return: min number of candidates
        :rtype: `int`

    .. automethod:: getCntCandidatesToKeepMax

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.accept_max`

        :return: max number of candidates
        :rtype: `int`

    .. automethod:: getCntMorphs

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.far_produce`

        :return:
        :rtype: `int`

    .. automethod:: getCntMorphsInDepth

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.close_produce`

        :return: number of morphs to generate
        :rtype: `int`

    .. automethod:: getDistToTargetDepthSwitch

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.far_close_threshold`

        :return: distance threshold
        :rtype: `float`

    .. automethod:: getCntMaxMorphs

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.max_morphs_total`

        :return: maximum number of non-producing morphs per molecule
        :rtype: `int`

    .. automethod:: getItThreshold

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.non_producing_survive`

        :return: number of generations until non-producing morphs are pruned
        :rtype: `int`

    .. automethod:: setSourceMol

        Set the `source molecule`.

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.source`

        :param \*args: `MolpherMol` instance or SMILES representing the current `source molecule`
        :type \*args: `MolpherMol` or `str`

    .. automethod:: setTargetMol

        Set the `target molecule`.

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.target`

        :param \*args: `MolpherMol` instance or SMILES representing the current `source molecule`
        :type \*args: `MolpherMol` or `str`

    .. automethod:: setChemOperators

        Set the current `chemical operators` as `str` identifiers.

        .. include:: oper_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.operators`

        :param \*args: current `chemical operators` as a `tuple` of `str` identifiers
        :type \*args: `tuple` of `str`

    .. automethod:: setFingerprint

        Set the `molecular fingerprint`.

        .. include:: fing_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.fingerprint`

        :param \*args: `molecular fingerprint` identifier
        :type \*args: `str`

    .. automethod:: setSimilarityCoef

        Set the `similarity measure`.

        .. include:: sim_table.rst

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.similarity`

        :param \*args: `similarity measure` identifier
        :type \*args: `str`

    .. automethod:: setMinAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.weight_min`

        :param \*args: minimum acceptable weight during filtering
        :type \*args: `float`

    .. automethod:: setMaxAcceptableMolecularWeight

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.weight_max`

        :param \*args: maximum acceptable weight during filtering
        :type \*args: `float`

    .. automethod:: setCntCandidatesToKeep

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.accept_min`

        :param \*args: min number of candidates
        :type \*args: `int`

    .. automethod:: setCntCandidatesToKeepMax

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.accept_max`

        :param \*args: max number of candidates
        :type \*args: `int`

    .. automethod:: setCntMorphs

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.far_produce`

        :param \*args:
        :type \*args: `int`

    .. automethod:: setCntMorphsInDepth

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.close_produce`

        :param \*args: number of morphs to generate
        :type \*args: `int`

    .. automethod:: setDistToTargetDepthSwitch

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.far_close_threshold`

        :param \*args: distance threshold
        :type \*args: `float`

    .. automethod:: setCntMaxMorphs

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.max_morphs_total`

        :param \*args: maximum number of non-producing morphs per molecule
        :type \*args: `int`

    .. automethod:: setItThreshold

        .. seealso:: `molpher.core.ExplorationParameters.ExplorationParameters.non_producing_survive`

        :param \*args: number of generations until non-producing morphs are pruned
        :type \*args: `int`