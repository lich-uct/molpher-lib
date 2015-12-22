Tutorial
========

..  todo:: write some very short intro (what the tutorial is going to be about and stuff)

Creating an Exploration Tree and Setting Morphing Parameters
------------------------------------------------------------

.. py:currentmodule:: molpher.core

As it was explained in the `../introduction`, Molpher maintains a data structure called an `exploration tree` to save and
evaluate putative `paths in chemical space <chemical space path>`. This data structure is represented by
a `molpher.core.ExplorationTree` instance in the Python API. Let's now initialize an exploration tree that will
represent a path between *cocaine* (a popular recreational drug) and *procaine* (compound that replaced cocaine
as a local anesthetic):

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :lines: 1,4-7
    :caption: The most basic way to initilize an exploration tree.
    :name: tree-init
    :linenos:

The code shown in  :numref:`tree-init` simply initilizes the tree from the supplied SMILES.
At the moment the tree is pretty simple. It only contains the root molecule (cocaine in this particular instance).
We can manipulate this instance and read data from it in multiple ways, but let's start by printing out the
`source molecule` and `target molecule` of the tree:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :lines: 9-10
    :caption: Printing out the source and target.
    :name: print-source-target
    :linenos:

If we run the code from :numref:`tree-init` and :numref:`print-source-target` in a script or on command line, it produces
the following output::

    loading SAScore.dat ... done
    Source:  CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    Target:  O=C(OCCN(CC)CC)c1ccc(N)cc1

just as we would expect.

.. note:: Besides our source and target, we can also see that a data file was loaded successfully. That means
    the `molpher` package was
    initilaized successfully and is ready for use.

The `ExplorationTree.params` dictionary doesn't just store the
source and target molecule, but also houses other `morphing parameters`. Let's take a look:

..  code-block:: python

    print(tree.params)

Output:

..  code-block:: python

    {
        'max_morphs_total': 1500,
        'far_close_threshold': 0.15,
        'weight_max': 100000.0,
        'close_produce': 150,
        'fingerprint': 'MORGAN',
        'accept_min': 50,
        'source': 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2',
        'target': 'O=C(OCCN(CC)CC)c1ccc(N)cc1',
        'weight_min': 0.0,
        'non_producing_survive': 5,
        'accept_max': 100,
        'operators': (
            'ADD_ATOM',
            'ADD_BOND',
            'BOND_CONTRACTION',
            'BOND_REROUTE',
            'INTERLAY_ATOM',
            'MUTATE_ATOM',
            'REMOVE_ATOM',
            'REMOVE_BOND'
        ),
        'far_produce': 80,
        'similarity': 'TANIMOTO'
    }

As we can see there are many more. We will explain the most important ones in this tutorial, but you can see the
description of `ExplorationParameters` (especially :numref:`param-table`) for a detailed reference.


We can adjust the morphing parameters during runtime as we like. All we need to do is just supply the `params` attribute
of our tree instance with a new dictionary:

..  code-block:: python

    tree.params = {
        'non_producing_survive' : 2
        'weight_max' : 500.0
    }

Output:

..  code-block:: python
    :emphasize-lines: 11, 4
    :linenos:

    {
        'max_morphs_total': 1500,
        'far_close_threshold': 0.15,
        'weight_max': 500.0,
        'close_produce': 150,
        'fingerprint': 'MORGAN',
        'accept_min': 50,
        'source': 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2',
        'target': 'O=C(OCCN(CC)CC)c1ccc(N)cc1',
        'weight_min': 0.0,
        'non_producing_survive': 2,
        'accept_max': 100,
        'operators': (
            'ADD_ATOM',
            'ADD_BOND',
            'BOND_CONTRACTION',
            'BOND_REROUTE',
            'INTERLAY_ATOM',
            'MUTATE_ATOM',
            'REMOVE_ATOM',
            'REMOVE_BOND'
        ),
        'far_produce': 80,
        'similarity': 'TANIMOTO'
    }

Here we just tighthened the constraints on molecular weight and decreased the number of acceptable 'non-producing'
`morph generations <morph generation>` to 2 (see :numref:`param-table` for details on the  morphing parameters).
If we supply an incomplete set of parameters (like in the above example), only the parameters specified in the
given dictionary will be changed.

..  warning:: Note that changing some parameters during runtime may have some adverse effects on the exploration
    process so use it with caution.

..  seealso:: The class `ExplorationParameters` can be used to hold parameters independently of a tree.
        It also includes a method that validates the parameters -- see `ExplorationParameters`
        and `ExplorationTree` for details on how to use this class.

Building an Exploration Tree
----------------------------

Building an exploration tree is a multi-step iterative process. This part of the tutorial outlines all the steps
involved and introduces the user to the most basic concepts of the library and molecular morphing.

Generating and Manipulating Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we know how to initilize an exploration tree and how to set morphing parameters, let's take a look at how the
`chemical space` exploration works in practice.

Let's generate some `morphs <morph>` from the current leaves of the tree first:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Generating and exploring morphs.
    :lines: 21-25
    :name: generate-morphs
    :linenos:

Output:

..  code-block:: none

    (<molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9ec2c0f0> >,)
    CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    (<molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9ec2c0f0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3253c0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3253f0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325420> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325450> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325480> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3254b0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3254e0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325510> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325540> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325570> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3255a0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3255d0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325600> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325630> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325660> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325690> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3256c0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3256f0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325720> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325750> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325780> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3257b0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3257e0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325810> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325840> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325870> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3258a0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3258d0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325900> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325930> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325960> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325990> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3259c0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c3259f0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325a20> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325a50> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325a80> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325ab0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325ae0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325b10> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325b40> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325b70> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325ba0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325bd0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325c00> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325c30> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325c60> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325c90> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325cc0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325cf0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325d20> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325d50> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325d80> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325db0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325de0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325e10> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325e40> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325e70> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325ea0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325ed0> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325f00> >, <molpher.swig_wrappers.core.MolpherMol; proxy of <Swig Object of type 'MolpherMol *' at 0x7f1c9c325f30> >)
    63

The code in :numref:`generate-morphs` first tells the tree to return its current leaves.
As we only have one molecule in the tree (our cocaine `source molecule`), the `leaves` member
only contains one element. We can verify that it is indeed our cocaine by calling its `getSMILES()` method.

The `generateMorphs()` method tells the tree to generate some `morphs <morph>`
from the current leaves for us (how many depends on the current state of the tree
and parameters `far_produce`, `close_produce` and `far_close_threshold`).

..  note:: The number of generated morphs also depends on other factors as well. For example, some structures
        might not be parsed correctly and hence might be removed from the list after the morphs are generated.

We can access the created morphs from the `candidates`
member of the tree instance. It is a `tuple` of `MolpherMol` instances. `MolpherMol` instances serve
as proxy objects to the underlying C++ representation of molecules generated by Molpher.

These instances can be used to read and manipulate compounds present in a tree.
However, because of their C++ origin, these objects cannot be treated as regular python objects.
The current implementation does not inform Python when and if the molecule is destroyed
(removed from the tree or otherwise invalidated). Therefore, one
simple rule have to be followed: **Do not save references to tree-bound molecules.**

You can tell whether a `MolpherMol` is tree-bound or not by calling its `isBound()` method.
If it returns `True`, the molecule is bound to a tree and, therefore, it is not safe to save it anywhere for later use.
If you want to have the information contained within a `MolpherMol` object available even when
the molecule is destroyed, you can tell it to replicate itself by calling its `copy()` method.
The replicated instance is not bound
to any tree and its lifetime is only affected by Python runtime and its garabage collector. However,
unbinding the molecule from the tree also means that the actual data in the tree will not be modified,
if this instance is changed. The following code example ilustrates this behaviour:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Copying morphs.
    :lines: 29-42
    :name: copying-morphs
    :linenos:

Output:

..  code-block:: none

    True
    0.90625
    0.5

    False
    0.5
    0.7
    0.5

We can see that when modifying a bound instance the value in the tree is modified as well,
but when we change an unbound instance the value stays the same.

..  note:: For more information on the available methods see the `MolpherMol` documentation.

Sorting and Filtering Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes the order of molecules in `candidates` might matter to us. As of yet, the only way to sort the
generated morphs in `candidates` is by calling the `sortMorphs()` method or using the `SortMorphsOper` operation
(see :ref:`operations` for more on operations). This sorts the molecules in the order of increasing distance
from the target. Placing the closest molecules at the very top of the list. For example, this order has meaning
for the built-in probability fileter (see `PROBABILITY` for details).

When the list of candidates is populated and sorted, we need to choose the morphs that
will form the next `generation <morph generation>`. The code below ilustrates
how we can sort the morphs and do the subsequent filtering manually:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Sorting and manually filtering morphs.
    :lines: 46-61
    :name: filtering-morphs
    :linenos:

Output:

..  code-block:: none

    [('COC(=O)C1C2CCC(CC1OC(=O)C1C=CC=C1)N2C', 0.5), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C', 0.7068965517241379), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143)]
    (True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True)
    63
    (True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False)

In :numref:`filtering-morphs`, `candidates_mask` member can be easily changed by writing
a `list` or `tuple` of new values into it. Here we simply select the first three morphs as the new `morph generation`.
Note that the list of selected morphs is sorted from the molecule closest to the target to the farthest,
because we called the `sortMorphs()` method previously.

..  warning:: The new mask must be the same length as the `candidates` member. If this requirement
        is not satisified, an instance of `RuntimeError` is raised.

..  warning:: The mask should only be set after the morphs are sorted. If the mask is set and
        the order of morphs is changed, the mask will stay the same and will have to be updated
        to follow the new order.

..  seealso:: The library also implements a few built-in filters. You can use the
        `filterMorphs()` method to invoke them. See the method's documentation for more information
        on the available filtering options.

Extending and Pruning
~~~~~~~~~~~~~~~~~~~~~

When we have the morphs we want to attach to the tree selected, we can call `extend()`
to connect them to their respective parents and make them the new leaves:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Extending the tree.
    :lines: 65-76
    :name: extending-tree
    :linenos:

Output:

..  code-block:: none

    [('COC(=O)C1C2CCC(CC1OC(=O)C1C=CC=C1)N2C', 0.5), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C', 0.7068965517241379), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143)]
    1
    False

We can see that after extending the tree, our selected morphs (see :numref:`filtering-morphs`) have become the new leaves and that the
`morph generation` counter `generation_count` was set to one. We could now repeat this process
and have Molpher iteratively explore further in `chemical space`. We can also read the `path_found` member
at every iteration to check for the presence of the `target molecule` in the tree and terminate the process,
if a path is found.

Because a tree generated in this way grows exponentially, a pruning strategy is needed in order
to keep the number of explored putative paths to a minimum by discarding those that are not getting any
closer to the `target molecule`.

We call the molecule that have not generated any morphs closer to the target than itself a *non-producing molecule*.
The number of generations to wait before removing the descendents of a `non-producing molecule`
from the tree is given by the `non_producing_survive` parameter.

Tree pruning can be requested anytime by calling the `prune()` method. In our example, the method didn't prune
any paths, because the `non_producing_survive` parameter is set to 2 generations in this particular instance.

..  seealso:: In addition to the `non_producing_survive` parameter, there is the `max_morphs_total` parameter,
        which imposes a restriction on the maximum number of
        descendents of one `non-producing molecule`. If the number of descendents
        reaches this threshold, the molecule is removed along with the descendents.

..  _operations:

Tree Operations
---------------

We call every action that is performed on an `exploration tree` a *tree operation*.
This concept is represented in the library as the `TreeOperation` abstract class.
This class becomes useful, for example, when we run into a situation where we need to build
several exploration trees at once. In such case we might want to apply the same set of operations
to every tree in the set. Moreover, this abstraction allows us to implement our own tree operations
and reuse code easily. We can run any operation on a given tree simply by supplying it to the
`runOperation()` method of the `ExplorationTree` class. The example below shows how to implement
the same exploration algorithm as demonstrated in the preceeding sections:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Using operations to implement simple chemical space exploration.
    :lines: 2-3,80-121
    :emphasize-lines: 3,27
    :name: operations-example
    :linenos:

Output:

..  code-block:: none

    1
    False
    [('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.7704918032786885), ('CN1C2CCC1C1C(=O)OCC2C1OC(=O)C1=CC=CC=C1', 0.7833333333333333)]

..  note:: Because the morphing algorithm is not deterministic, the set of obtained morphs is slightly different from the
        one in the above examples.

In :numref:`operations-example` we see an example where we create a list of operations to perform on a tree. Most of
them are built-in operations (discussed below), but we also chose to define our own operation for the filtering step
(see the highlighted lines).
In order to do that, we simply create a subclass of the `TreeOperation` abstract class and we override its
:meth:`~molpher.swig_wrappers.core.TreeOperation.__call__()` method with the implementation we need --
even though that our new filter is very likely to cause the exploration to get stuck in local optima.

..  note:: In the above example we also redefine the `TreeOperation.getTree()` method so that it 'casts' the basic
        `molpher.swig_wrappers.core.ExplorationTree` (returned by `TreeOperation.getTree()`)
        to the 'more pythonic' `molpher.core.ExplorationTree`.
        This override will probably be integrated in the library itself soon. Here, it just serves as another example
        of how we can change the default behaviour.

Each operation can have a tree associated with it, but it is not necessary.
We can verify if a tree is associated with an operation by calling its :meth:`~molpher.swig_wrappers.core.TreeOperation.getTree()`
method. If there is no tree associated with the instance, it returns `None`.

..  note:: Most of the built-in operations will raise a `RuntimeError`, if invoked without a tree attached to them.

Built-in Operations
~~~~~~~~~~~~~~~~~~~

A few operations are already defined in the library:

    - `GenerateMorphsOper`
    - `SortMorphsOper`
    - `FilterMorphsOper`
    - `FindLeavesOper`
    - `ExtendTreeOper`
    - `PruneTreeOper`
    - `TraverseOper`

They are all dervied from `TreeOperation` and contain the full set of operations performed on a tree in
the original Molpher algorithm as published in [1]_. Therefore, the original algorithm can be
implemented using those operations.

In this tutorial, we will pay particular attention to the `TraverseOper` operation, because it differs
from the others and introduces the concept of a tree callback function (see `tree-traversal`).
For more details on each operation, see the particular pages in the documentation.

.. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

..  _tree-traversal:

Traversing the Tree
^^^^^^^^^^^^^^^^^^^

A special place among the operations belongs to the `TraverseOper` class. It does not directly implement a part
of a morphing algorithm, but serves as a means of traversing molecules in a tree and reading/modifying them
as needed. Let's ilustrate this with an example:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Traversing the tree using a callback.
    :lines: 125-138
    :name: traverse-example
    :linenos:

Output:

..  code-block:: none

    # Root #
    SMILES:  CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    Descendents:  ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=CC=C2)CCN1C', 'CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 'COC(=O)C1C2C(=O)CC(CC1OCC1=CC=CC=C1)N2C')
    # Morph #
    Parent: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    SMILES:  COC(=O)C1C2C(=O)CC(CC1OCC1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    SMILES:  CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    SMILES:  CCC1C(C(=O)OC)C(OC(=O)C2=CC=CC=C2)CCN1C
    Descendents:  ()

In :numref:`traverse-example`, we derive from the `TraverseCallback` class, which is an abstract class with
an abstract method called :meth:`~molpher.swig_wrappers.core.TraverseCallback.processMorph()`. This method takes one argument,
which is a `MolpherMol` instance
of a molecule in the tree. We need to override this method in our derived class in order to implement our own
behaviour.

The callback is then associated with a `TraverseOper` instance, which can be run on a tree as any other
tree operation. When the operation is run it traverses the tree from root to leaves and injects
every molecule it encounters into our :meth:`~molpher.swig_wrappers.core.TraverseCallback.processMorph()` method.

..  note:: We can also pass a molecule to the `TraverseOper` constructor. In that case, a subtree will be traversed
        using the specified molecule as the root of the subtree.

There is also a much more convenient way to traverse the tree. The `ExplorationTree` class implements
the `traverse()` method. It simply takes any python callable and tries to use it instad of the
:meth:`~molpher.swig_wrappers.core.TraverseCallback.processMorph()` method.
However, under the hood it does the same thing as we did in :numref:`traverse-example`.
Therefore, the above code can be turned into:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Traversing the tree using a callback -- the simple version.
    :lines: 142-151
    :name: short-traverse-example
    :linenos:

Output:

..  code-block:: none

    # Root #
    SMILES:  CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    Descendents:  ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=CC=C2)CCN1C', 'CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 'COC(=O)C1C2C(=O)CC(CC1OCC1=CC=CC=C1)N2C')
    # Morph #
    Parent: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    SMILES:  COC(=O)C1C2C(=O)CC(CC1OCC1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    SMILES:  CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    SMILES:  CCC1C(C(=O)OC)C(OC(=O)C2=CC=CC=C2)CCN1C
    Descendents:  ()

Summary
-------

..  todo:: write

