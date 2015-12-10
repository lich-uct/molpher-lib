TreeOperation
=============

.. py:currentmodule:: molpher.swig_wrappers.core

.. autoclass:: TreeOperation

    This class defines an interface for `tree operations <tree operation>`.
    It can be executed on any `exploration tree` using the `ExplorationTree.runOperation()` method
    or a tree can be explicitly attached to it.

    It is an abstract class and therefore cannot be instantiated, but it defines a constructor that
    either takes a single argument (an instance of `ExplorationTree`) or nothing at all (in which case
    the associated `exploration tree` instance is set to `None`).

    .. automethod:: __call__

        This method is abstract and is meant to be everriden by derived classes to
        run the operation on the currently registered tree (if any).

    ..  automethod:: getTree

        Fetch the `ExplorationTree` instance attached to this operation.

        .. note:: This method *always* returns an instance of `molpher.swig_wrappers.core.ExplorationTree`.
            If you want to change this behaviour, you will need to override this method in the derived class.

        :return: an `exploration tree`
        :rtype: `ExplorationTree`

    ..  automethod:: setTree

        Associate the given tree with this operation instance.

        :param tree: the tree to associate this instance with
        :type tree: `ExplorationTree`