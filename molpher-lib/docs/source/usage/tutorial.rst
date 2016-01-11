Tutorial
========

This section gives an overview of most of the features currently available in the library
in the form of a very simple example. The example script is located in the `examples` package
and can be downloaded from :download:`here <../../../python/molpher/examples/basics.py>`. The used example
`XML template` is available from :download:`this link <../../../python/molpher/examples/cocaine-procaine-template.xml>`
(see `templates-snapshots` for more details on templates). If you want to know more about the library or
the source code, we suggest you check out the `source-code-docs`.

.. contents:: Table of Contents

Creating an Exploration Tree and Setting Morphing Parameters
------------------------------------------------------------

.. py:currentmodule:: molpher.core

As it was explained in the `../introduction`, Molpher maintains a data structure called an `exploration tree` to save and
evaluate putative `paths in chemical space <chemical space path>`. This data structure is represented by
a `molpher.core.ExplorationTree` instance in the Python API. This tutorial will show how to create and build
an exploration tree that will try to search for a path
between *cocaine* (a popular recreational drug) and *procaine* (compound that replaced cocaine
as a local anesthetic). Both of these compounds act on sodium channels in the neuronal cell membrane
and probably have the same mode of action, which makes them ideal candidates for a morphing experiment.
Let's now initialize our exploration tree:

..  code-block:: python
    :caption: The most basic way to initialize an exploration tree.
    :name: tree-init
    :linenos:

    from molpher.core import ExplorationTree as ETree

    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    tree = ETree(source=cocaine, target=procaine)

The code shown in  :numref:`tree-init` simply initializes the tree from the supplied SMILES.
At the moment the tree is pretty simple. It only contains the root molecule (cocaine in this particular instance).
We can manipulate this instance and read data from it in multiple ways, but let's start by printing out the
`source molecule` and `target molecule` of the tree:

..  code-block:: python
    :caption: Printing out the source and target.
    :name: print-source-target
    :linenos:

    print('Source: ', tree.params['source'])
    print('Target: ', tree.params['target'])

If we run the code from :numref:`tree-init` and :numref:`print-source-target` in a script or on command line, it produces
the following output::

    loading SAScore.dat ... done
    Source:  CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    Target:  O=C(OCCN(CC)CC)c1ccc(N)cc1

.. note:: Besides the information about our source and target, we can also see that a data file was loaded successfully. That means
    the :mod:`molpher` package was
    initialized successfully and is ready for use.

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
    print(tree.params)

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

Here we just tightened the constraints on molecular weight and decreased the number of acceptable 'non-producing'
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

Now that we know how to initialize an exploration tree and how to set morphing parameters, let's take a look at how the
`chemical space` exploration works in practice.

Let's generate some `morphs <morph>` from the current leaves of the tree first:

..  code-block:: python
    :caption: Generating and reading morphs.
    :name: generate-morphs
    :linenos:

    print(tree.leaves)
    print(tree.leaves[0].getSMILES())
    tree.generateMorphs()
    print(tree.candidates)
    print(len(tree.candidates))

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
to any tree and its lifetime is only affected by Python runtime and its garbage collector. However,
unbinding the molecule from the tree also means that the actual data in the tree will not be modified,
if this instance is changed. The following code example illustrates this behaviour:

..  code-block:: python
    :caption: Copying morphs.
    :name: copying-morphs
    :linenos:

    candidate = tree.candidates[0]
    print(candidate.isBound())
    print(tree.candidates[0].getDistToTarget())
    candidate.setDistToTarget(0.5)
    print(tree.candidates[0].getDistToTarget())

    print()

    candidate_copy = candidate.copy()
    print(candidate_copy.isBound())
    print(candidate_copy.getDistToTarget())
    candidate_copy.setDistToTarget(0.7)
    print(candidate_copy.getDistToTarget())
    print(candidate.getDistToTarget())

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
for the built-in probability filter (see `PROBABILITY` for details).

When the list of candidates is populated and sorted, we need to choose the morphs that
will form the next `generation <morph generation>`. The code below illustrates
how we can sort the morphs and do the subsequent filtering manually:

..  code-block:: python
    :caption: Sorting and manually filtering morphs.
    :name: filtering-morphs
    :linenos:

    tree.sortMorphs()
    print(tree.candidates_mask)
    print(len(tree.candidates_mask))
    mask = [False for x in tree.candidates_mask]
    mask[0] = True
    mask[1] = True
    mask[2] = True
    tree.candidates_mask = mask
    print(tree.candidates_mask)
    print(
        [
            (x.getSMILES(), x.getDistToTarget())
            for idx,x in enumerate(tree.candidates)
            if tree.candidates_mask[idx]
        ]
    )

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
        is not satisfied, an instance of `RuntimeError` is raised.

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

..  code-block:: python
    :caption: Extending the tree.
    :name: extending-tree
    :linenos:

    tree.extend()
    print(
        sorted(
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )
    print(tree.generation_count)
    print(tree.path_found)
    tree.prune()

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
This concept is represented in the library as the :class:`molpher.core.operations.TreeOperation` abstract class.
This class becomes useful, for example, when we run into a situation where we need to build
several exploration trees at once. In such case we might want to apply the same set of operations
to every tree in the set. Moreover, this abstraction allows us to implement our own tree operations
and reuse code easily. We can run any operation on a given tree simply by supplying it to the
`runOperation()` method of the `ExplorationTree` class. The example below shows how to implement
the same exploration algorithm as demonstrated in the preceding sections:

..  code-block:: python
    :caption: Using operations to implement simple chemical space exploration.
    :emphasize-lines: 3,17
    :name: operations-example
    :linenos:

    from molpher.core.operations import *

    class MyFilterMorphs(TreeOperation):

        def __call__(self):
            mask = [False for x in self.tree.candidates_mask]
            mask[0] = True
            mask[1] = True
            mask[2] = True
            self.tree.candidates_mask = mask

    tree = ETree(source=cocaine, target=procaine)

    iteration = [
        GenerateMorphsOper()
        , SortMorphsOper()
        , MyFilterMorphs()
        , ExtendTreeOper()
        , PruneTreeOper()
    ]

    for oper in iteration:
        tree.runOperation(oper)

    print(tree.generation_count)
    print(tree.path_found)
    print(
        sorted(
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )

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
In order to do that, we simply create a subclass of the :class:`~molpher.core.operations.TreeOperation` abstract class and we override its
:meth:`~molpher.swig_wrappers.core.TreeOperation.__call__()` method with the implementation we need --
even though that our new filter is very likely to cause the exploration to get stuck in local optima.

Each operation can have a tree associated with it, but it is not necessary.
We can verify if a tree is associated with an operation by calling its :meth:`~molpher.swig_wrappers.core.TreeOperation.getTree()`
method or accessing the `TreeOperation.tree` attribute of the class. If there is no tree associated with the instance, they both return `None`.

..  note:: The built-in operations will raise a `RuntimeError`, if invoked without a tree attached to them.

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

They are all dervied from :class:`~molpher.core.operations.TreeOperation` and contain the full set of operations performed on a tree in
the original Molpher algorithm as published in [1]_. Therefore, the original algorithm can be
implemented using those operations.

In this tutorial, we will pay particular attention to the `TraverseOper` operation, because it differs
from the others and introduces the concept of a tree callback function (see `tree-traversal`).
For more details on each operation, see the designated pages in the documentation.

.. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

..  _tree-traversal:

Traversing the Tree
^^^^^^^^^^^^^^^^^^^

A special place among the operations belongs to the `TraverseOper` class. It does not directly implement a part
of a morphing algorithm, but serves as a means of traversing molecules in a tree and reading/modifying them
as needed. Let's illustrate this with an example:

..  code-block:: python
    :caption: Traversing the tree using a callback.
    :name: traverse-example
    :linenos:

    class MyCallback(TraverseCallback):

        def processMorph(self, morph):
            if not morph.getParentSMILES():
                print("# Root #")
            else:
                print('# Morph #')
                print('Parent:', morph.getParentSMILES())
            print('SMILES: ', morph.getSMILES())
            print('Descendents: ', morph.getDescendants())

    callback = MyCallback()
    traverse = TraverseOper(callback)
    tree.runOperation(traverse)

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
the `traverse()` method. It simply takes any python callable and tries to use it instead of the
:meth:`~molpher.swig_wrappers.core.TraverseCallback.processMorph()` method.
However, under the hood it does the same thing as we did in :numref:`traverse-example`.
Therefore, the above code can be turned into:

..  code-block:: python
    :caption: Traversing the tree using a callback -- the simple version.
    :name: short-traverse-example
    :linenos:

    def process(morph):
        if not morph.getParentSMILES():
            print("# Root #")
        else:
            print('# Morph #')
            print('Parent:', morph.getParentSMILES())
        print('SMILES: ', morph.getSMILES())
        print('Descendents: ', morph.getDescendants())

    tree.traverse(process)

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

..  _templates-snapshots:

Templates and Tree Snapshots
----------------------------

We don't always have to initialize `morphing parameters` by hand. We can use a `XML template` instead.
Here is an example of a template file:

..  literalinclude:: ../../../python/molpher/examples/cocaine-procaine-template.xml
    :language: xml
    :caption: A complete XML template file.
    :name: template-file
    :linenos:

A `XML template` is similar to a configuration file and can be loaded like an ordinary snapshot (see :numref:`loading-snapshot`).
The following example shows loading of a `XML template`, creating a tree from it, extending the tree and saving
a tree snapshot:

..  code-block:: python
    :caption: Loading a template and saving a built tree as a XML snapshot.
    :name: saving-snapshot
    :linenos:

    template_file = 'cocaine-procaine-template.xml'

    tree = ETree.createFromSnapshot(template_file)
    print(tree.params)

    for oper in iteration:
        tree.runOperation(oper)

    print(
        sorted(
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )

    tree.saveSnapshot('snapshot.xml')

Output:

..  code-block:: none

    Parse molecule CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2 >> COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Parse molecule O=C(OCCN(CC)CC)c1ccc(N)cc1 >> CCN(CC)CCOC(=O)C1=CC=C(N)C=C1
    The new iteration has been created from template:
    source: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    target: CCN(CC)CCOC(=O)C1=CC=C(N)C=C1

    Snapshot successfully created from: cocaine-procaine-template.xml
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
    [('COC(=O)C(COC(=O)C1=CC=CC=C1)C1CCC(C)N1C', 0.7627118644067796), ('CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.7704918032786885), ('COC(=O)C1C2CCC(CC1COC(=O)C1=CC=CC=C1)N2C', 0.7741935483870968)]

We can see that all the parameters are the same as in the `XML template` and that
the resulting tree can be built using the same list of operations
as in :numref:`operations-example`. We even get the same set of newly generated leaves.
In :numref:`saving-snapshot` we also want to serialize our tree instance to disk so we save it as
a snapshot using the `saveSnapshot()` method.

The saved tree can be later reconstructed with the
:meth:`~molpher.core.ExplorationTree.ExplorationTree.createFromSnapshot()` factory method:

..  warning:: Note that there is currently a bug inside Molpher that prevents
        loading of descendents of molecules in the tree. Therefore, the code
        from :numref:`loading-snapshot` does not work correctly at the moment.

..  code-block:: python
    :caption: Loading a snapshot of a previously generated tree.
    :name: loading-snapshot
    :linenos:

    new_tree = ETree.createFromSnapshot('snapshot.xml')
    print(new_tree.params)
    print(
        sorted(
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in new_tree.leaves
        ], key=lambda x : x[1]
        )
    )

Output:

..  code-block:: none

    Snapshot successfully created from: snapshot.xml
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
    [('COC(=O)C(COC(=O)C1=CC=CC=C1)C1CCC(C)N1C', 0.7627118644067796), ('CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.7704918032786885), ('COC(=O)C1C2CCC(CC1COC(=O)C1=CC=CC=C1)N2C', 0.7741935483870968)]

Example Exploration Algorithm Implementations
---------------------------------------------

Let's now wrap up this tutorial with two example implementations of a morphing experiment (see :numref:`complete-example`
and :numref:`bidirectional-example`).

Using the bits of code we have created above. The script in :numref:`complete-example` shows how to implement
a search for a path in `chemical space` between *cocaine* and *procaine* with a customized filtering step:

..  literalinclude:: ../../../python/molpher/examples/experiment.py
    :language: python
    :caption: Example implementation of a pathfinding algorithm.
    :name: complete-example
    :linenos:

Output:

..  code-block:: none

    Loading data from: /home/sichom/Projects/Molpher/Molpher_repo/molpher-lib/python/molpher/swig_wrappers/SAScore.dat
    loading SAScore.dat ... done
    Parse molecule CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2 >> COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Parse molecule O=C(OCCN(CC)CC)c1ccc(N)cc1 >> CCN(CC)CCOC(=O)C1=CC=C(N)C=C1
    Iteration 1
    [('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C', 0.7068965517241379), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Iteration 2
    [('COC(=O)CC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5769230769230769), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=CC=C1', 0.6415094339622642), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Iteration 3
    [('CCC(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)OC', 0.4897959183673469), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=C(N)C=C1', 0.4897959183673469), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=CC=C1', 0.6415094339622642), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Iteration 4
    [('CCC(CCOC(=O)C1=CC=C(N)C=C1)N(CC)CCC(=O)O', 0.44680851063829785), ('CCC(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.4565217391304348), ('CCC(CC(=O)O)N(C)CCCOC(=O)C1=CC=C(N)C=C1', 0.46808510638297873), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=CC=C1', 0.6415094339622642), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Iteration 5
    [('CCN(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.3571428571428571), ('CCC(CC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.40909090909090906), ('CCN(CCC(=O)O)C(C)CCOC(=O)C1=CC=C(N)C=C1', 0.4347826086956522), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=CC=C1', 0.6415094339622642), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Pruned (descendents only): COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1
    Iteration 6
    [('CN(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.3157894736842105), ('CCN(CCCOC(=O)C1=CC=C(N)C=C1)CCC(=O)O', 0.31707317073170727), ('CCN(CCC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.3414634146341463), ('CCC(CC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.40909090909090906), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=CC=C1', 0.6415094339622642), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Candidate morph: CCN(CCCOC(=O)C1=CC=C(N)C=C1)CCC(=O)O already present in the tree. Skipping...
    Pruned (descendents only): COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C
    Pruned (descendents only): CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C
    Iteration 7
    [('CCN(CCCO)CCCOC(=O)C1=CC=C(N)C=C1', 0.30000000000000004), ('CN(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.3157894736842105), ('CCN(CCCOC(=O)C1=CC=C(N)C=C1)CC(=O)O', 0.31707317073170727), ('CCN(CCC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.3414634146341463), ('CCC(CC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.40909090909090906), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Pruned (descendents only): COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C
    Iteration 8
    [('CCN(CCOC(=O)C1=CC=C(N)C=C1)CC(=O)O', 0.18918918918918914), ('CCN(CCCO)CCOC(=O)C1=CC=C(N)C=C1', 0.21052631578947367), ('CCN(CCCOC(=O)C1=CC=C(N)C=C1)CCOO', 0.275), ('CN(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.3157894736842105), ('CCN(CCC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.3414634146341463), ('CCC(CC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.40909090909090906), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Iteration 9
    [('CCCN(CC)CCOC(=O)C1=CC=C(N)C=C1', 0.1428571428571429), ('CCN(CCOO)CCOC(=O)C1=CC=C(N)C=C1', 0.16666666666666663), ('CCN(CCO)CCOC(=O)C1=CC=C(N)C=C1', 0.16666666666666663), ('CCN(CCOC(=O)C1=CC=C(N)C=C1)CC(=O)O', 0.18918918918918914), ('CN(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.3157894736842105), ('CCN(CCC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.3414634146341463), ('CCC(CC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.40909090909090906), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Candidate morph: CCCN(CC)CCOC(=O)C1=CC=C(N)C=C1 already present in the tree. Skipping...
    Pruned (descendents only): CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=C(N)C=C1
    Iteration 10
    [('CCN(CC)CCOC(=O)C1=CC=C(N)C=C1', 0.0), ('CCN(CO)CCOC(=O)C1=CC=C(N)C=C1', 0.1428571428571429), ('CCN(CCOO)CCOC(=O)C1=CC=C(N)C=C1', 0.16666666666666663), ('CCN(CCOC(=O)C1=CC=C(N)C=C1)CC(=O)O', 0.18918918918918914), ('CN(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)O', 0.3157894736842105), ('CCN(CCC(=O)O)N(C)CCOC(=O)C1=CC=C(N)C=C1', 0.3414634146341463), ('CCC(CC(=O)OC)N(C)CCCOC(=O)C1=CC=C(N)C=C1', 0.4897959183673469), ('COCCC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C', 0.5416666666666667), ('CCC1C(C(=O)OC)C(OC(=O)C2=CC=C(N)C=C2)CCN1C', 0.6666666666666667), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('COC(=O)C1C2CC1N(C)C2CCOC(=O)C1=CC=CC=C1', 0.728813559322034)]
    Path found:
    COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C 1.7976931348623157e+308
    COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C 0.7068965517241379
    COC(=O)CC1CCC(CCOC(=O)C2=CC=C(N)C=C2)N1C 0.5769230769230769
    CCC(CCOC(=O)C1=CC=C(N)C=C1)N(C)CCC(=O)OC 0.4897959183673469
    CCC(CCOC(=O)C1=CC=C(N)C=C1)N(CC)CCC(=O)O 0.44680851063829785
    CCN(CCC(=O)O)C(C)CCOC(=O)C1=CC=C(N)C=C1 0.4347826086956522
    CCN(CCCOC(=O)C1=CC=C(N)C=C1)CCC(=O)O 0.31707317073170727
    CCN(CCCO)CCCOC(=O)C1=CC=C(N)C=C1 0.30000000000000004
    CCN(CCCO)CCOC(=O)C1=CC=C(N)C=C1 0.21052631578947367
    CCCN(CC)CCOC(=O)C1=CC=C(N)C=C1 0.1428571428571429
    CCN(CC)CCOC(=O)C1=CC=C(N)C=C1 0.0

    Process finished with exit code 0

The above implementation is nothing more than just the tutorial code bits inside a loop. The loop checks if a path was found at each iteration.
If the path is found it is backtracked through the tree and printed out as a sequence of molecules from source to target.

The second example we have here is a little more elaborate, but implements a very simple idea. Instead of one exploration tree,
we build two trees that each searches for a path to the closest molecule in the other tree:

..  literalinclude:: ../../../python/molpher/examples/bidirectional.py
    :language: python
    :caption: Example implementation of a bidirectional pathfinding algorithm.
    :name: bidirectional-example
    :linenos:

Output (only the final print of the path is shown):

..  code-block:: python

    [
        'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C',
        'CCC1C(C(=O)OC)C(OC(=O)C2=CC=CC=C2)CCN1C',
        'CCC1C(COC)C(OC(=O)C2=CC=CC=C2)CCN1C',
        'CCC1C(COC)C(OC(=O)C2=CC=CC=C2)CCN1CC',
        'CCN1CCC(OC(=O)C2=CC=CC=C2)C(COC)C1C',
        'CCN1CC(OC(=O)C2=CC=CC=C2)C(COC)C1C',
        'CCN1CC(OC(=O)C2=CC=CC=C2)C(CO)C1C',
        'CCC1C(C)N(CC)CC1OC(=O)C1=CC=CC=C1',
        'CCC1C(C)C(OC(=O)C2=CC=CC=C2)CN1CC',
        'CCC1CC(OC(=O)C2=CC=CC=C2)CN1CC',
        'CCC1C(OC(=O)C2=CC=CC=C2)CN1CC',
        'CCC1C(OC(=O)C2=CC=C(N)C=C2)CN1CC',
        'CCN1CC(OC(=O)C2=CC=C(N)C=C2)C1C',
        'CCN(CC)CCOC(=O)C1=CC=C(N)C=C1'
    ]
    Total Execution Time: 153964.77300000002 # in milliseconds

This 'bidirectional search' algorithm uses the built-in operations to facilitate the search,
but does one extra procedure after
an iteration is completed -- it changes the target molecules of the trees.

When the new leaves are connected, both trees are traversed and molecules
closest to the current target are identified in each. The closest molecule from one tree is then
set as the new target for the tree searching in the opposite direction and vice versa.

In :numref:`bidirectional-example` we also use the `time.clock` function to measure the execution
times of each potentially time-consuming operation.

Summary
-------

If you have read the tutorial all the way to here, you now probably have a decent idea on what this library does
and how to use it. If you have any suggestions on how to improve the library or bug reports, please send them to
sichom@vscht.cz. All help on the project is much appreciated.
For more information on the source code itself refer to the `source-code-docs`.

