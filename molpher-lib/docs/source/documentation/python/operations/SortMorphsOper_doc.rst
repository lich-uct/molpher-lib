SortMorphsOper
==============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: SortMorphsOper
    :show-inheritance:

    A `tree operation` used to sort the `candidate morphs` in an `exploration tree`
    according to their distance to the target molecule at the moment of their creation.

    `Candidate morphs <candidate morphs>` are sorted in order of their increasing distance from the target.

    .. automethod:: __call__

        Sort the `candidate morphs` according to their distance (asceding order)
        to the target molecule at the moment of their creation.

