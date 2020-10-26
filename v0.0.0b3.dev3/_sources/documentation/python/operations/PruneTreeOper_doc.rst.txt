PruneTreeOper
=============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: PruneTreeOper
    :show-inheritance:

    This :term:`tree operation` prunes some unwanted molecules from the tree. Unwanted molecules
    are those that haven't generated any :term:`morphs <morph>` closer to the target in a given number
    of iterations (see `ExplorationData.non_producing_survive`) or those that produced way too many morphs
    without any of them being closer to the target (see `ExplorationData.max_morphs_total`).

    .. automethod:: __call__

        Prune the attached `ExplorationTree` instance according to the rules given in `ExplorationData`.

