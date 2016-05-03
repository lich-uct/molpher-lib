ExtendTreeOper
==============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: ExtendTreeOper
    :show-inheritance:

    A :term:`tree operation` that commits a :term:`morphing iteration` of the associated :term:`exploration tree` by attaching
    the :term:`candidate morphs` that survived the filtering to the tree so that they can become the new
    leaves and can be used to create the next generation of morphs.

    .. automethod:: __call__

        Attach the :term:`candidate morphs` in `ExplorationTree.candidates`
        (not masked in `ExplorationTree.candidates_mask`) as the new leaves of the tree.

