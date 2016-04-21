FilterMorphsOper
================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: FilterMorphsOper
    :show-inheritance:

    :param \*args: constructor takes either no arguments at all or the following positional arguments (in that order):
    :param tree: the tree to operate on (optional)
    :param filters: the union of filtering options to use (optional, see the descriptions below for more information)
    :param verbose: toggles if a verbose output is required (optional, defaults to `False`)

    This operation can filter the `candidate morphs` in an `exploration tree`
    according to a given set of filtering options that it implements. The current set of
    available filters is in the table below. For more detailed information on each filter
    see the member descriptions of this class.

    ..  seealso:: `ExplorationTree.filterMorphs()` and `molpher.core.ExplorationData`

    ..  include:: ../filters_table.rst

    .. autoattribute:: PROBABILITY

        If the number of derived morphs for a molecule exceeds `ExplorationData.accept_min`,
        then the probability of accepting every extra molecule is determined
        according to the following formula:
        :math:`\frac{0.25 - (idx - morphs_{min})}{4(morphs_{all} - morphs_{min})}`, where :math:`morphs_{min}`
        is the `ExplorationData.accept_min` parameter, :math:`morphs_{all}` is the current
        number of molecules in the `exploration tree` and :math:`idx` is the position of the given `morph`
        in `ExplorationTree.candidates`.

    .. autoattribute:: WEIGHT

        Uses the `ExplorationData.weight_min` and `ExplorationData.weight_max` options
        to remove `candidate morphs` that do not satisfy the weight constraints.

    .. autoattribute:: SYNTHESIS

        .. todo:: say something about the synthetic feasibility in more detail

    .. autoattribute:: MAX_DERIVATIONS

        Each `ExplorationTree` instance keeps track of the number of morphs generated from each
        molecule in the tree. If there are `candidate morphs` in `ExplorationTree.candidates`
        such that more than `ExplorationData.max_morphs_total` number of morphs would be derived from
        the given molecule, then these `candidate morphs` will be filtered out.

    .. autoattribute:: DUPLICATES

        Remove all duplicate morphs (same canonical SMILES representation
        as a molecule already present in the tree).

        .. note:: This filter is always turned on and can't be disabled at the moment.

    .. autoattribute:: HISTORIC_DESCENDENTS

        The tree also keeps track of all molecules previously derived from any molecule in the tree.
        This filter removes `candidate morphs` that have alredy been previously derived from their parent molecule.

    .. autoattribute:: ALL

        Use all of the filters mentioned above.
