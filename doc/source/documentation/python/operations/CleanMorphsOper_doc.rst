CleanMorphsOper
===============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: CleanMorphsOper
    :show-inheritance:

    :param \*args: the :term:`exploration tree` to clean
    :type \*args: `ExplorationTree` (optional)

    This operation removes all morphs from `candidates` that
    are marked for removal in `candidates_mask`.

    .. automethod:: __call__

        Run this operation on the attached tree.

