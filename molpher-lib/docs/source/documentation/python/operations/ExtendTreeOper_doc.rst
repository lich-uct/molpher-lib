ExtendTreeOper
==============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExtendTreeOper
    :show-inheritance:

    A `tree operation` that commits a `morphing iteration` of the associated `exploration tree` by attaching
    the `candidate morphs` that survived the filtering to the tree so that they can become the new
    leaves and can be used to create the next generation of morphs.

    .. automethod:: __call__

        Attach the `candidate morphs` in `ExplorationTree.candidates`
        (not masked in `ExplorationTree.candidates_mask`) as the new leaves of the tree.

