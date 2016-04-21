FindLeavesOper
==============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: FindLeavesOper
    :show-inheritance:

    A `tree operation` used to access the leaves of an `exploration tree`.

    .. automethod:: __call__

        Fetch the leaves from the attached tree.

    .. automethod:: fetchLeaves

        Returns the leaves of the last tree on which `__call__` method was applied.

        :return: the leaves of the attached tree
        :rtype: `tuple` of :py:class:`~molpher.core.MolpherMol.MolpherMol` instances