MolpherMol
==========

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: MolpherMol

    :param other: another instance of :py:class:`MolpherMol` or a SMILES string
    :type other: :py:class:`MolpherMol` or `str`

    .. note:: If another :py:class:`MolpherMol` instance is specified, the
        tree ownership is not preserved in the new instance.
        It just acts as a place holder for the data and these are not updated if
        the original instance changes.

    Represents a molecule (usually attached to an :term:`exploration tree`). It allows
    the user of the library to read data about the molecule as well as modify it.

    .. automethod:: copy

        Creates a new :py:class:`MolpherMol` instance, which is not bound to any tree.
        It has the same effect as calling the constructor with another instance
        as its argument.

        .. attention:: Any modifications to the returned molecule will
            have no effect on the original instance.

        :return:  new instance
        :rtype: :py:class:`MolpherMol`

    .. automethod:: getSMILES

        Returns canonical SMILES of the molecule.

        :return:  canonical SMILES representation of the molecule
        :rtype: :py:class:`str`

    .. automethod:: getDistToTarget

        Returns the distance of this molecule to the target molecule.

        By default computed as :math:`1 - similarity`, where
        :math:`similarity` is the similarity computed using the specified :term:`molecular fingerprint`
        and :term:`similarity measure` (see `molpher.core.selectors` for
        more information on the available :term:`selectors`).

        .. note:: Distance to the target can be easily modified, if needed (see `setDistToTarget`).

        .. seealso:: `molpher.core.ExplorationData`

        :return:  distance to target molecule
        :rtype: :py:class:`float`

    .. automethod:: getParentSMILES

        Canonical SMILES of the direct parent molecule in the respective :term:`exploration tree`.

        :return:  canonical SMILES of parent
        :rtype: :py:class:`str`

    .. automethod:: getDescendants

        Returns all direct descendants of this molecule in the respective :term:`exploration tree`
        as canonical SMILES. Therefore, all molecules returned by this method
        should be present in the tree at the time of the call.

        :return:  an iterable of descendants
        :rtype: :py:class:`tuple` of :py:class:`str` objects

    .. automethod:: getHistoricDescendants

        Returns the historic descendants (all morphs created from this molecule)
        of this molecule in the respective :term:`exploration tree` as canonical SMILES.
        Some of the molecules might have been removed from the tree before the
        call so some of them might be missing from the tree already.

        :return:  an iterable of historic descendants
        :rtype: :py:class:`tuple` of :py:class:`str` objects

    .. automethod:: getItersWithoutDistImprovement

        Returns the number of morph generations derived from this molecule
        that did not contain any morphs with an improvement in structural
        distance from the :term:`target molecule`.

        .. seealso:: :attr:`~molpher.core.ExplorationData.ExplorationData.non_producing_survive`,
                :attr:`~molpher.core.MolpherMol.MolpherMol.gens_without_improvement`

        :return:  number of generations without distance improvement
        :rtype: :py:class:`int`

    .. automethod:: getSAScore

        The synthesis feasibility score.

        .. todo:: specify more accurately

        :return:  synthesis feasibility score
        :rtype: :py:class:`float`

    .. automethod:: getMolecularWeight

        Returns the molecular weight of this compound.

        :return:  molecular weight
        :rtype: :py:class:`float`

    .. automethod:: setSAScore

        Set the synthesis feasibility score
        as computed using :term:`SAScore.dat`.

        :param score:  new synthesis feasibility score
        :type score: :py:class:`float`

    .. automethod:: setDistToTarget

        Set the distance to the :term:`target molecule`.

        :param dist:  new distance to target
        :type dist: :py:class:`float`

    .. automethod:: getTree

        Return the :term:`exploration tree` associated with
        this molecule.

        :param dist:  new distance to target
        :type dist: :py:class:`ExplorationTree` or `None`

    .. automethod:: setItersWithoutDistImprovement

        Sets the number of morph generations derived from this molecule
        that did not contain any morphs with an improvement in structural
        distance from the :term:`target molecule`.

        .. seealso:: :attr:`~molpher.core.ExplorationData.ExplorationData.non_producing_survive`,
                :attr:`~molpher.core.MolpherMol.MolpherMol.gens_without_improvement`
