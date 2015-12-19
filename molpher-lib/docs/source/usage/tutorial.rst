Tutorial
========

..  todo:: write some very short intro (what the tutorial is going to be about and stuff)

Creating an Exploration Tree
----------------------------

.. py:currentmodule:: molpher.core

As it was explained in the `../introduction`, Molpher maintains a data structure called an `exploration tree` to save and
evaluate putative `paths in chemical space <chemical space path>`. This data structure is represented by
a `molpher.core.ExplorationTree` instance in the Python API. Let's now initialize an exploration tree that will
represent a path between *cocaine* (a popular recreational drug) and *procaine* (compound that replaced cocaine
as a local anesthetic):

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :lines: 1-6
    :caption: The most basic way to initilize an exploration tree.
    :name: tree-init
    :linenos:

The code shown in  :numref:`tree-init` simply initilizes the tree from the supplied SMILES.
At the moment the tree is pretty simple. It only contains the root molecule (cocaine in this particular instance).
We can manipulate this instance and read data from it in multiple ways, but let's start by printing out the
`source molecule` and `target molecule` of the tree:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :lines: 8-9
    :caption: Printing out the source and target.
    :name: print-source-target
    :linenos:

If we run the code from :numref:`tree-init` and :numref:`print-source-target` in a script or command line, it produces
the following output::

    loading SAScore.dat ... done
    Source:  CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    Target:  O=C(OCCN(CC)CC)c1ccc(N)cc1

This will print out the source and target molecule in the way that we would expect.

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
of our instance with a new dictionary:

..  code-block:: python

    tree.params = {
        'non_producing_survive' : 2,4
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

Extending the Exploration Tree
------------------------------

Extending an exploration tree is a multi-step process. This part of the tutorial outlines all the steps
involved and introduces the user to the most basic concepts of the library.

Generating and Manipulating Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we know how to initilize an exploration tree and how to set morphing parameters, let's take a look at how the
`chemical space` exploration works in practice.

Let's generate some `morphs <morph>` from the current leaves of the tree first:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Generating and exploring morphs.
    :lines: 20-24
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

These instances can be used to read and manipulate compounds present in the tree.
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
unbinding the molecule from the tree also means that the actual data in the tree will not be modified
if this instance is changed. The following code example ilustrates this behaviour:

..  literalinclude:: ../../../python/molpher/examples/basics.py
    :language: python
    :caption: Copying morphs.
    :lines: 28-41
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

Filtering Morphs
~~~~~~~~~~~~~~~~

..  todo:: write

Tree Operations
---------------

..  todo:: about tree operations