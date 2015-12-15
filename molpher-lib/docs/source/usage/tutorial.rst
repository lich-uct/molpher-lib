Tutorial
========

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

Here we just tighthened the constraints on molecular weight and decrease the number of acceptable 'non-producing'
`morph generations <morph generation>` to 2 (see :numref:`param-table` for details on the  morphing parameters).
If we supply an incomplete set of parameters (like in the above example), only the parameters specified in the
given dictionary will be changed.

..  warning:: Note that changing some parameters during runtime may have some adverse effect on the exploration
    process so use this approach with caution.