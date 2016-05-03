GenerateMorphsOper
==================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: GenerateMorphsOper
    :show-inheritance:

    A :term:`tree operation` used to generate the :term:`candidate morphs` in an :term:`exploration tree`.

    .. automethod:: __call__

        Generate new :term:`candidate morphs` and save them to the `candidates` member of
        the attached `ExplorationTree` instance.

