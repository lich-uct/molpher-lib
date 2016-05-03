SortMorphsOper
==============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: SortMorphsOper
    :show-inheritance:

    A :term:`tree operation` used to sort the :term:`candidate morphs` in an :term:`exploration tree`
    according to their distance to the target molecule at the moment of their creation.

    :term:`Candidate morphs <candidate morphs>` are sorted in order of their increasing distance from the target.

    .. automethod:: __call__

        Sort the :term:`candidate morphs` according to their distance (ascending order)
        to the target molecule at the moment of their creation.

