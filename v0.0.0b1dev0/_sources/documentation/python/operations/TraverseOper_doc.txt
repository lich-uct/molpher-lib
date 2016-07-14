TraverseOper
============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: TraverseOper
    :show-inheritance:

    :param \*args: the :term:`exploration tree` and the callback to call. Optionally we can specify
        a :term:`morph` which will define the root of the subtree to traverse. If no root
        is specified, then the :term:`source molecule` is selected automatically.
    :type \*args: `ExplorationTree`, `TraverseCallback` and :py:class:`~molpher.core.MolpherMol.MolpherMol` (optional)

    A :term:`tree operation` that can register a `TraverseCallback` and call its
    :py:meth:`~TraverseCallback.__call__` method on every molecule it encounters in the given subtree
    (defined by its root).

    ..  warning:: The callback to Python is realized using `SWIG's director feature
            <http://www.swig.org/Doc3.0/Python.html#Python_directors>`_,
            which makes it possible to keep the implementation
            concurrent and efficient. However, there is a problem:

            If a reference to the `MolpherMol`
            instance is saved into a variable that outlives the call to the
            callback function, this reference then becomes invalid when the call is
            finished. Therefore, doing something like this:

            .. code-block:: python

                var = None
                class MyCallback(TraverseCallback):

                    def __call__(self, morph):
                        if var:
                            print("Previous:", var.smiles)
                        var = morph

                callback = MyCallback()
                traverse = TraverseOper(callback)
                tree.runOperation(traverse)

            will likely result in a segmentation fault upon a second call to the callback function,
            because the object referenced by ``var`` will no longer refer to valid memory.
            This is likely a result of SWIG freeing the pointer without taking
            the existing reference from Python into account.

            One way to work around this is
            to just use the SMILES string to fetch the
            reference directly using the :meth:`~molpher.swig_wrappers.core.ExplorationTree.fetchMol` method:

            .. code-block:: python

                var = morph.tree.fetchMol(morph.smiles)

    .. automethod:: __call__

        Run the attached callback on every molecule in the subtree with the given root.

