TraverseCallback
================

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: TraverseCallback

    A simple callback that can process morphs in an `exploration tree`.
    When supplied to `TraverseOper` the morphs that belong to the currently traversed
    subtree are processed using this method.

    .. note:: This is an abstract class. You will have to declare a derived class to use its features.

    .. seealso:: `TraverseOper`

    .. automethod:: processMorph

       An abstract method that should be overriden in derived classes.
       It should implement the morph processing.

       :param morph: currently processed morph
       :type morph: `MolpherMol`