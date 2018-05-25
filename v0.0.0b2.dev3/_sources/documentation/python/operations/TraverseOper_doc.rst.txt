TraverseOper
============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: TraverseOper
    :show-inheritance:

    :param \*args: the :term:`exploration tree` and the callback to call. Optionally, a canonical SMILES string
        of a :term:`morph` which will define the root of the subtree to traverse can be specified. If no root
        is given, then the :term:`source molecule` is selected automatically and the traversal will be performed on the whole tree.
    :type \*args: `ExplorationTree`, `TraverseCallback` and `str` (optional)

    A :term:`tree operation` that can register a `TraverseCallback` and call its
    :py:meth:`~TraverseCallback.__call__` method on every molecule it encounters in the given subtree
    (defined by its root).

    .. automethod:: __call__

        Run the attached callback on every molecule in a subtree with the given root.

