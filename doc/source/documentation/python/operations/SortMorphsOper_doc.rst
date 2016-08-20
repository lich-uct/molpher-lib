SortMorphsOper
==============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: SortMorphsOper
    :show-inheritance:

    :param tree: a tree to run this operation on (optional)
    :type tree: :class:`~molpher.swig_wrappers.core.ExplorationTree` instance
    :param callback: callable to call every time two molecules need to be compared
    :type callback: instance of :class:`~molpher.swig_wrappers.core.SortMorphsCallback`

    A :term:`tree operation` used to sort the :term:`candidate morphs` in an :term:`exploration tree`
    using the specified :class:`~molpher.swig_wrappers.core.SortMorphsCallback`.

    :term:`Candidate morphs <candidate morphs>` are sorted in order of their increasing distance from the target.

    .. automethod:: __call__

        Sort the :term:`candidate morphs`.

