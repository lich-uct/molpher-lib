ExplorationTree
===============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExplorationTree

    :param \*args: one positional argument (either an `molpher.swig_wrappers.core.ExplorationParameters`
        instance or a `str` (SMILES representing the `source molecule`)
    :type \*args: `str` or `molpher.swig_wrappers.core.ExplorationParameters`

    Represents an `exploration tree` and facilitates the interaction with it.

    .. automethod:: createFromSnapshot(snapshot)

        Returns an `ExplorationTree` instance using a supplied `ExplorationTreeSnapshot`.

        :param snapshot: an `ExplorationTreeSnapshot` to generate an `ExplorationTree` from
        :type snapshot: `ExplorationTreeSnapshot`
        :return: `ExplorationTree` instance created from the supplied snapshot
        :rtype: `ExplorationTree`

    .. automethod:: createSnapshot

        Creates an `ExplorationTreeSnapshot` instance that can be saved to disk and contains information about
        the molecules currently present in the `exploration tree` along with the morphing parameters.

        .. note:: Only the molecules attached to the cuurent `exploration tree` are serialized.
            If this instance contains unattached candidate morphs, this data will be lost.

        :return: `ExplorationTreeSnapshot` instance created from this tree
        :rtype: `ExplorationTreeSnapshot`

    .. automethod:: runOperation

        Takes a `TreeOperation` instance and attempts to run it by attaching this tree
        instance to it and calling the `TreeOperation.__call__()` method.

        .. seealso:: `operations/classes`

        .. note:: Doing this will also automatically attach the given `TreeOperation` instance to
            this `ExplorationTree` object.

        :param operation: a `TreeOperation` instance to run
        :type operation: `TreeOperation`

    .. automethod:: fetchLeaves

        Gives a `tuple` of `MolpherMol` instances that represent the leaves of this `ExplorationTree`.

        :return: `tuple` of `MolpherMol` instances representing the leaves of this tree
        :rtype: `tuple`

    .. automethod:: fetchMol(canonSMILES)

        Attempts to fetch a molecule (represented by its canonical SMILES string) from the current tree.
        If the molecule exists, it returns
        an appropriate `MolpherMol` instance. Otherewise, it raises a `RuntimeError`.

        :param canonSMILES: a `str` representing SMILES of the molecule being searched for
        :type canonSMILES: `str`
        :return: `MolpherMol` representing the looked up molecule
        :rtype: `MolpherMol`
        :raises: `RuntimeError` if the molecule specified could not be found in the tree

    .. automethod:: hasMol(canonSMILES)

        Returns `True` or `False` to indicate the presence of a molecule in the current tree.

        :param canonSMILES: a `str` representing SMILES of the molecule being searched for
        :type canonSMILES: `str`
        :return: `bool` value representing the presence of the searched molecule in this tree (`True` if it is present, `False` otherwise)
        :rtype: `bool`

    .. automethod:: deleteSubtree

        Deletes a specified molecule from the tree and all of its descendents.

        :param canonSMILES: a `str` representing SMILES of the molecule being removed
        :type canonSMILES: `str`

    .. automethod:: generateMorphs

        Generates mew `candidate morphs` using the current `exploration parameters`
        and saves them (see `ExplorationParameters` for details).

        The generated compounds can be retrieved using the `getCandidateMorphs` method.

    .. automethod:: sortMorphs

        Sorts the `candidate morphs` according to their distances from the `target molecule`.

        ..  note:: When generated each `morph` is assigned a distance to current target.
                This distance *is not* updated if the target of the tree changes in the future.

    .. automethod:: filterMorphs

        Filters the `candidate morphs` according to the supplied filtering option (see `FilterMorphsOper`
        for a list of available filters). The options can be easily combined and passed
        to the method using the ``|`` operator. For example, the method can be called like this:
        ``tree.filterMorphs(FilterMorphsOper.SYNTHESIS | FilterMorphsOper.PROBABILITY)``. If no
        filter information is passed, all filters are used automatically
        (equal to passing `FilterMorphsOper.ALL`).

        The results of the filtering can be observered using the `getCandidateMorphsMask` method.

        :param \*args: one or more filtering options combined using the ``|`` operator
        :type \*args: `FilterMorphsOper` options

    .. automethod:: extend

        Attach all `candidate morphs` that have not been filtered out to the tree.
        In other words, make them the new leaves.

    .. automethod:: prune

        Perform the pruning of the tree according to the rules specified by current
        `exploration parameters` (see `ExplorationParameters` for details).

    .. automethod:: getThreadCount

        Returns the number of threads this particular `ExplorationTree` will use. The value of 0
        indicates that the same number of threads as the number of available cores will be used.

    .. automethod:: getParams

        Can be used to obtain the current `exploration parameters`.

        :return: instance of `ExplorationParameters`, which holds the current configuration options for the computations on the tree
        :rtype: `ExplorationParameters`

    .. automethod:: getCandidateMorphs

        Gives a `tuple` of `MolpherMol` instances that represent the `candidate morphs`
        of this `ExplorationTree` instance.

        :return: `tuple` of `MolpherMol` instances representing the `candidate morphs`
        :rtype: `tuple`

    .. automethod:: getCandidateMorphsMask

        Returns a mask for the `tuple` returned by `getCandidateMorphs`. The filtered out morphs
        are denoted by `False` at the respective position in the mask.

        Morphs filtered out in this way will not be attached to the tree when `extend` is called.

        :return: `tuple` of `bool` instances which shows what `candidate morphs` were filtered out
        :rtype: `tuple`

    .. automethod:: setThreadCount

        Sets the maximum number of threads used by this instance. The value of 0
        indicates that the same number of threads as the number of available cores will be used.

        :param threadCnt: an `int` indicating the number of threads to use
        :type threadCnt: `int`

    .. automethod:: setParams

        Change the current `exploration parameters` for this instance.

        ..  warning:: This may invalidate some data in the tree (such as the distances from the `target molecule` computed so far).
                Use with caution.

        :param params: an instance of `ExplorationParameters`, which represents the new `exploration parameters` for this instance
        :type params: `int`

    .. automethod:: setCandidateMorphsMask(mask)

        Sets the morphs' mask of this instance to a custom value. If the length of the new mask does not
        correspond to the number of `candidate morphs` in the tree, a `RuntimeError` will be raised.

        :param mask: `tuple` of `bool` instances which shows what `candidate morphs` will be filtered out
            (not attached to the tree when `extend` is called)
        :type mask: `tuple`
        :raises: `RuntimeError` if the length of the new mask does not correspond to the number of `candidate morphs` in the tree