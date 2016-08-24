SortMorphsCallback
==================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: SortMorphsCallback

    A simple callback that can be used to override how morphs are sorted.
    When supplied to `SortMorphsOper`, the `__call__` method is executed
    every time two morphs need to be compared by the sorting algorithm.
    The return value of `__call__` indicates the result of the comparison.

    .. note:: This is an abstract class. You will have to declare a derived class to that implements
        the abstract features.

    .. seealso:: `SortMorphsOper`

    .. automethod:: __call__

       An abstract method overridden in derived classes.
       It should return `True` or `False`. `True` indicates that the
       :samp:`morph_a` is greater than the :samp:`morph_b`.

       :param morph_a: first morph in comparison
       :type morph_a: :py:class:`MolpherMol`
       :param morph_b: second morph in comparison
       :type morph_b: :py:class:`MolpherMol`
       :rtype: `bool`