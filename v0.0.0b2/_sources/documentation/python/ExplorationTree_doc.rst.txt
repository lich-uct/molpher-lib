ExplorationTree
===============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExplorationTree

    Represents an `exploration tree` and facilitates the most basic interaction with it. It
    is a tight wrapper around the C++ implementation. The `molpher.core.ExplorationTree` class
    extends upon this implementation and provides a more user-friendly version of this class.

    ..  attention:: This class does not define a constructor. Use
            the `create()` factory method to spawn instances.

    .. automethod:: create(*args)

        Returns an `ExplorationTree` instance using a supplied `ExplorationData`,
        :term:`tree snapshot` or SMILES strings
        of the :term:`source molecule` and the :term:`target molecule`.

        :param \*args: an `ExplorationData` instance, path to a :term:`tree snapshot`
            or a :term:`source molecule` and a :term:`target molecule` as two arguments.
        :type \*args: `ExplorationData` or `str`
        :return: `ExplorationTree` instance created from the supplied snapshot
        :rtype: `ExplorationTree`

    .. automethod:: save(filename)

        Saves the tree to a snapshot file. The filename should
        end with either ``.xml`` or ``.snp``.

        :param filename: path to the new snapshot file
        :type filename: `str`

    .. automethod:: runOperation

        Takes a `TreeOperation` instance and attempts to run it by attaching itself
        to it and calling its :py:meth:`~TreeOperation.__call__()` method.

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

        Generates mew :term:`candidate morphs` using the current :term:`exploration parameters`
        and saves them (see `ExplorationData` for details).

        The generated compounds can be retrieved using the `getCandidateMorphs` method.

        ..  note:: All morphs already present in the tree are filtered out automatically.

    .. automethod:: sortMorphs

        Sorts the :term:`candidate morphs` according to their distances from the :term:`target molecule`.

        ..  warning:: When morphs are generated the distance to the current target is saved for each morph.
                This distance *is not* updated if the target of the tree changes afterwards.

    .. automethod:: filterMorphs

        Filters the :term:`candidate morphs` according to the supplied filtering options
        (see the table below for a listing of the available filters).
        The filtering options can be easily combined and passed
        to the method using the ``|`` operator. For example, the method can be called like this:
        ``tree.filterMorphs(FilterMorphsOper.SYNTHESIS | FilterMorphsOper.PROBABILITY)``. If no
        filter information is passed, all filters are used automatically
        (equal to passing `FilterMorphsOper.ALL`).

        ..  include:: filters_table.rst

        The results of the filtering can be observed using the `getCandidateMorphsMask` method.

        ..  seealso:: `FilterMorphsOper`

        :param \*args: one or more filtering options combined using the ``|`` operator
        :type \*args: `FilterMorphsOper` options

    .. automethod:: extend

        Attach all :term:`candidate morphs` that have not been filtered out to the tree.
        In other words, make them the new leaves.

        ..  note:: By calling this method a `morphing iteration` is committed.

    .. automethod:: prune

        Perform the pruning of the tree according to the rules specified by current
        :term:`exploration parameters` (see `ExplorationData` for details).

    .. automethod:: getThreadCount

        Returns the number of threads this particular `ExplorationTree` will use. The value of 0
        indicates that the same number of threads as the number of available cores will be used.

    .. automethod:: asData

        Can be used to obtain the current :term:`exploration parameters`.

        :return: instance of `ExplorationData`, which holds the current configuration options for the computations on the tree
        :rtype: `ExplorationData`

    .. automethod:: getCandidateMorphs

        Gives a `tuple` of `MolpherMol` instances that represent the `candidate morphs`
        of this `ExplorationTree` instance.

        :return: `tuple` of `MolpherMol` instances representing the `candidate morphs`
        :rtype: `tuple`

    .. automethod:: getCandidateMorphsMask

        Returns a mask for the `tuple` returned by `getCandidateMorphs`. The filtered out morphs
        are denoted by `False` at the respective position in the mask.

        Morphs filtered out in this way will not be attached to the tree when :meth:`extend` is called.

        :return: `tuple` of `bool` instances which shows what :term:`candidate morphs` were filtered out
        :rtype: `tuple`

    .. automethod:: getGenerationCount

        Returns the number of :term:`morph generations <morph generation>` in the tree. Basically
        the number of :term:`morphing iterations <morphing iteration>` performed so far.

        :return: number of :term:`morphing iterations <morphing iteration>` performed so far
        :rtype: `int`

    .. automethod:: setThreadCount

        Sets the maximum number of threads used by this instance. The value of 0
        indicates that the same number of threads as the number of available cores will be used.

        :param threadCnt: an `int` indicating the number of threads to use
        :type threadCnt: `int`

    .. automethod:: update(data)

        Change the current :term:`exploration parameters` for this instance. If
        the supplied `ExplorationData` instance contains a tree topology, it is
        ignored and only the parameters are taken into account.

        ..  warning:: This may invalidate some data in the tree (such as the distances from
                the :term:`target molecule` computed so far).

        :param param: an instance of `ExplorationData`, which represents the new :term:`exploration parameters` for this instance
        :type param: `ExplorationData`

    .. automethod:: setCandidateMorphsMask(mask)

        Sets the morphs' mask of this instance to a custom value. If the length of the new mask does not
        correspond to the number of :term:`candidate morphs` in the tree, a `RuntimeError` will be raised.

        :param mask: `tuple` of `bool` instances which shows what :term:`candidate morphs` will be filtered out
            (not attached to the tree when `extend()` is called)
        :type mask: `tuple`
        :raises: `RuntimeError` if the length of the new mask does not correspond to
            the number of :term:`candidate morphs` in the tree

    .. automethod:: traverse(*args)

        Traverses the tree and calls the ``()`` of the callback
        on each morph in it. A SMILES string of a root
        of a subtree can be specified before the callback in order
        to traverse just the subtree.

        :param \*args: either an instance of a class derived from `TraverseCallback` or
                both SMILES of the root of a subtree and a `TraverseCallback`
        :type \*args: `TraverseCallback` or `str` and `TraverseCallback`