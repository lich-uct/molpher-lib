Tutorial
========

This section gives a comprehensive overview of both the
morphing algorithm itself as well as the
software features currently available
By writing and commenting on an example Python script (located in the `examples` package
and available for download :download:`here <../../../src/python/molpher/examples/basics.py>`)
we present the building blocks that anyone can easily put together to implement their own :term:`Molpher`.

However, not all of the features and settings will be presented and we encourage
you to take a look at the `source-code-docs` if you want to know more about the implementation itself
or more advanced topics.

..  note:: We do not provide usage examples or a description of the C++ API in this tutorial,
        but the Python interface actually follows the C++ implementation very closely
        and the classes from the `core` package are in fact derived from the proxy
        classes generated automatically by `SWIG <http://www.swig.org/>`_ when the C++ API is
        wrapped for Python.

.. contents:: Table of Contents

.. _morphing-algorithm:

The Morphing Algorithm
----------------------

Let us first shortly describe how :term:`molecular morphing` works. This should help anyone who has never heard
about the method to understand it and get an idea on what this library can do for them.

In order to find a path that connects the :term:`source molecule` and the :term:`target molecule`
(see `../introduction` if you are not familiar with these terms),
Molpher maintains a data structure called an :term:`exploration tree` to save and
evaluate possible :term:`paths in chemical space <chemical space path>`. It is basically a
*directed rooted tree* (with the edges directed away from the root) where the source molecule acts as the root vertex.

During :term:`molecular morphing`, the tree is grown by modifying the current leaves
with predefined :term:`chemical operators`. These operators act on
an input structure by producing a new molecule as output, which is structurally very close to its parent compound.
We can also see chemical operators as the labels of the edges of the tree that
connect two compounds between which a transformation occurred, the parent and the child.

By connecting the newly created molecules to the tree as the children of their parents, we can
iteratively generate multiple :term:`chemical space paths <chemical space path>` at once and evaluate
them according to a certain objective function. At the moment the objective function is
simply the structural distance between the newest molecule on a path and the target molecule, but
the user of the library is free to define his own criteria. For example, he/she can easily implement
a method that will only extend paths that lead to compounds satisfying a given pharmacophore or
physicochemical properties.

By continuously morphing the compounds in this way, we can effectively 'travel' through :term:`chemical space`
towards various areas of interest. Thanks to the flexible API of the library this 'journey' can be realized
in many ways and can accommodate almost any exploration strategy one might think of. For example,
in this tutorial we will, among other things,
show how to :ref:`implement our own building blocks to use with the library <operations>` or
:ref:`make two trees cooperate with each other as they search for a single path <bidirectional>`.

Creating an Exploration Tree and Setting Morphing Parameters
------------------------------------------------------------

.. py:currentmodule:: molpher.core

Using the library the tree is represented by
an instance of the `molpher.core.ExplorationTree` class and in this tutorial
we will use it to search for a :term:`chemical space path`
between *cocaine* (a famous recreational drug) and *procaine* (compound that replaced cocaine
as local anesthetic). Both of these compounds act as blockers of the sodium channels in the neuronal cell membrane
which results in their anesthetic effects.

Let us now explore the :term:`chemical space` 'between' cocaine and procaine
by combining their structural features using the *Molpher-lib* library.
We will initialize the exploration tree first:

..  code-block:: python
    :caption: The most basic way to initialize an :term:`exploration tree` in Python.
    :name: tree-init
    :linenos:

    from molpher.core import ExplorationTree as ETree

    cocaine = 'CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2'
    procaine = 'O=C(OCCN(CC)CC)c1ccc(N)cc1'

    tree = ETree.create(source=cocaine, target=procaine) # initialize a tree that searches for a path from cocaine to procaine

The code shown in  :numref:`tree-init` simply initializes the tree from the supplied SMILES.
At the moment the tree is pretty simple. It only contains the root molecule (cocaine in this particular instance).
We can manipulate this instance and read data from it in multiple ways, but let's start by printing out the
:term:`source molecule` and :term:`target molecule` of the tree:

..  code-block:: python
    :caption: Printing out the :term:`source molecule` and :term:`target molecule` of a tree.
    :name: print-source-target
    :linenos:

    # print the smiles of the source and target
    print('Source: ', tree.params['source'])
    print('Target: ', tree.params['target'])

If we run the code from both :numref:`tree-init` and :numref:`print-source-target` in a script or command line, it produces
the following output:

..  code-block:: none

    loading SAScore.dat ... done
    Source:  COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Target:  CCN(CC)CCOC(=O)C1=CC=C(N)C=C1

Notice that the SMILES strings we supplied to the :py:meth:`~molpher.core.ExplorationTree.ExplorationTree.create`
method of the :py:class:`~molpher.core.ExplorationTree.ExplorationTree` class are not the same as those
that we supplied. That is because *Molpher-lib* automatically generates canonical SMILES for every
molecule that is created.

..  attention:: The library also discards any information about the stereocenters in the molecules,
    because the current implementation does not account for stereochemistry and treats all enantiomers
    as the same molecule. You might want to this in mind when working with the library.

..  note:: Besides the information about our source and target, we can also see that a data file was loaded successfully.
    That means the :mod:`molpher` package was initialized successfully and is ready for use. The data file itself is
    used to compute the synthetic feasibility scores for the generated morphs
    (can be read from any molecule by reading its `sascore` attribute).

The `ExplorationTree.params` dictionary does not just store the
source and target, but also houses other `morphing parameters`. Let's take a look:

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
        'source': 'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C',
        'target': 'CCN(CC)CCOC(=O)C1=CC=C(N)C=C1',
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

As we can see there is quite a lot of parameters that we can set.
We will explain the most important ones in this tutorial, but you can see the
documentation for the :py:class:`~ExplorationData.ExplorationData` class
(especially :numref:`param-table`) for a detailed reference.

We can adjust the morphing parameters during runtime as we like. All we need to do is just supply the `params` attribute
of our tree instance with a new dictionary:

..  code-block:: python

    # change selected parameters using a dictionary
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
        'source': 'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C',
        'target': 'CCN(CC)CCOC(=O)C1=CC=C(N)C=C1',
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
:term:`morph generations <morph generation>` to 2 (see :numref:`param-table` to get more information on what this parameter
does). If we supply an incomplete set of parameters (like in the above example),
only the parameters specified in the given dictionary will be changed.

..  warning:: Changing individual values in the `params` dictionary will have no effect.
    You always need to store a dictionary instance in it. This is because the value
    is regenerated every time the attribute is accessed to always reflect the current
    set of parameters valid for the current instance.

..  seealso:: :py:class:`~ExplorationData.ExplorationData`

Generating Morphs and Extending the Exploration Tree
----------------------------------------------------

This part of the tutorial outlines all the steps
involved in generating new morphs from the current leaves of the tree.
It also describes how the generated morphs can be filtered using the value of
the default objective function (structural distance from the target molecule),
the tree extended and the unfavorable paths (or their parts) removed from the tree.

Generating and Manipulating Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we know how to initialize an exploration tree and how to set morphing parameters,
we will take a look at how the :term:`chemical space` exploration works in practice.

Let us generate a few `morphs <morph>` from the current leaves of the tree first:

..  code-block:: python
    :caption: Generating and reading morphs.
    :name: generate-morphs
    :linenos:

    print(tree.leaves) # show the current leaves of the tree (only the source so far)
    print(tree.leaves[0].getSMILES())
    tree.generateMorphs() # generate new morphs
    print(tree.candidates)
    print(len(tree.candidates))

Output:

..  code-block:: none

    (<molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bc30>,)
    COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Generated 62 morphs.
    (<molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bc30>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bc00>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bc90>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bcc0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bcf0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bd20>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bd50>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bd80>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bdb0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bde0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78be10>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78be40>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78be70>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bea0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bed0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bf00>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bf30>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bf60>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bf90>, <molpher.core.MolpherMol.MolpherMol at 0x7f22ea78bfc0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06030>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06060>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06090>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae060c0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae060f0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06120>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06150>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06180>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae061b0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae061e0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06210>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06240>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06270>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae062a0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae062d0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06300>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06330>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06360>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06390>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae063c0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae063f0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06420>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06450>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06480>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae064b0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae064e0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06510>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06540>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06570>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae065a0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae065d0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06600>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06630>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06660>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06690>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae066c0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae066f0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06720>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06750>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae06780>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae067b0>, <molpher.core.MolpherMol.MolpherMol at 0x7f22eae067e0>)
    62

The code in :numref:`generate-morphs` first tells the tree to return its current leaves.
As we only have one molecule in the tree (our cocaine :term:`source molecule`), the `leaves` member
only contains one element. We can verify that it is indeed our cocaine by asking
the underlying :py:class:`~MolpherMol.MolpherMol` instance for its SMILES
using the `smiles` attribute.

..  hint:: The :py:class:`~MolpherMol.MolpherMol` class has a lot of useful attributes that
    can often be written into as well. You can easily replace the computed value of the objective
    function with your own, for example. See the documentation for the :py:class:`~MolpherMol.MolpherMol`
    class to get an overview of how you can modify the molecules in the tree or the generated candidate morphs.

The `generateMorphs()` method tells the tree to generate some `morphs <morph>`
from the current leaves for us. How many morphs will be generated depends
mostly on the current state of the tree
and parameters `far_produce`, `close_produce` and `far_close_threshold`.
However, it also depends on other factors. For example, some structures
might not be parsed correctly and hence might be removed from the list after the morphs are generated.
Also, a different number of morphs can be generated each time the method is run. That si due to
the non-deterministic character of the morphing algorithm which chooses the morphing operators to
use and parts of the structure to modify randomly.

We can access the newly generated morphs from the `candidates`
member of the tree instance. It is a `tuple` of :py:class:`~MolpherMol.MolpherMol` instances.
:py:class:`~MolpherMol.MolpherMol` instances serve
as proxy objects for the underlying C++ representation of molecules generated by :term:`Molpher`.
Hence, these instances can be used to read and manipulate compounds currently present in a tree
or the generated morphs.

..  attention:: The molecules saved in the `candidates` attribute of the tree actually do not
        belong to the tree just yet. See :ref:`extend-prune` for more information on
        how it is with tree ownership of molecules.

You can tell any molecule to replicate itself by calling its :py:meth:`~MolpherMol.MolpherMol.copy()` method.
The replicated instance is never bound to any tree and all the changes made to the copy are only made to the copy
and not the original instance. This is useful when we want to save the molecule from a tree for later use
and make sure that its state does not change as we keep growing the tree in our program.

The following code example illustrates how we can change and copy :py:class:`~MolpherMol.MolpherMol` instances:

..  code-block:: python
    :caption: Modifying and copying candidate morphs.
    :name: copying-morphs
    :linenos:

    # get the first morph in the candidate list
    candidate = tree.candidates[0]
    # print distance to target
    print(tree.candidates[0].dist_to_target)
    # set new distance to target
    candidate.dist_to_target = 0.5
    # look in the list of candidates and print new distance
    print(tree.candidates[0].dist_to_target)

    print()

    # make a copy of our molecule
    candidate_copy = candidate.copy()
    # set a new distance for the copy and verify that the original was not affected
    print(candidate_copy.dist_to_target)
    candidate_copy.dist_to_target = 0.7
    print(candidate_copy.dist_to_target)
    print(candidate.dist_to_target)
    print(tree.candidates[0].dist_to_target)

Output:

..  code-block:: none

    0.828125
    0.5

    0.5
    0.7
    0.5
    0.5

Sorting and Filtering Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes the order of the newly generated molecules in the `candidates` list might have a meaning to us;
for example, this order is used by the built-in probability filter (see `PROBABILITY` for details)
to compute the probabilities of survival for each candidate.

As of yet, the only way to sort the
generated morphs is by calling the :py:meth:`~molpher.swig_wrappers.core.ExplorationTree.sortMorphs`
method on the tree instance or using the :py:class:`~operations.SortMorphsOper.SortMorphsOper` operation
(see :ref:`operations` for more). This sorts the molecules in the order of increasing value
of the objective function (distance from the :term:`target molecule` by default).

Let us now sort the candidate morphs in our tree:

..  code-block:: python
    :caption: Sorting morphs according to the value of the objective function.
    :linenos:

    # sort the candidates in the tree according to their distance from target
    tree.sortMorphs()

    # verify
    print(tree.candidates_mask)
    print(
        [
            (x.smiles, x.dist_to_target)
            for idx,x in enumerate(tree.candidates)
        ]
    )

Output:

..  code-block:: none

    (True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True)
    [('COC(=O)C1C2CCC(CC1OC(=O)C1=CCC(C)C=C1)N2C', 0.5), ('CCC1N(C)C2CC(OC(=O)C3=CC=CC=C3)C12C(=O)OC', 0.7868852459016393), ('COC(=O)C1C2C3CN2C(C3)CC1OC(=O)C1=CC=CC=C1', 0.7868852459016393), ('CCC1CC(OC(=O)C2=CC=CC=C2)C(C(=O)OC)CN1C', 0.7868852459016393), ('COC(=O)C1C2C(=O)CC(CC1OCC1=CC=CC=C1)N2C', 0.7903225806451613), ('CCN1C2CCC1C(C(=O)OC)C(OC(=O)C1=CC=CC=C1)C2', 0.7936507936507937), ('COC(=O)C1C2CCC(C(N)C1OC(=O)C1=CC=CC=C1)N2C', 0.8032786885245902), ('COCC1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8064516129032258), ('COC(=O)C1C2CCC(C1OC(=O)C1=CC=CC=C1)N2C', 0.8103448275862069), ('COC(=O)CC(CC1CCCN1C)OC(=O)C1=CC=CC=C1', 0.8125), ('COC(=O)C1C2CCC(CNC1OC(=O)C1=CC=CC=C1)N2C', 0.8125), ('COC(=O)C1C2CCC(CC1(N)OC(=O)C1=CC=CC=C1)N2C', 0.8125), ('CONC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8125), ('CNC1CCCCC(OC(=O)C2=CC=CC=C2)C1C(=O)OC', 0.8135593220338984), ('COC(=O)C1C2CC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8166666666666667), ('COC(=O)C1C(C)N(C)C(C)CC1OC(=O)C1=CC=CC=C1', 0.8166666666666667), ('COC(=O)C1C2CCC(N2)C(C)C1OC(=O)C1=CC=CC=C1', 0.8166666666666667), ('CN1C2CCC1C(C(=O)O)C(OC(=O)C1=CC=CC=C1)C2', 0.8166666666666667), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)C2C', 0.819672131147541), ('COC(=O)C1(OC(=O)C2=CC=CC=C2)CC2CCC1N2C', 0.819672131147541), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2', 0.819672131147541), ('CC1CC2CC(OC(=O)C3=CC=CC=C3)C(C(=O)O)C1N2C', 0.819672131147541), ('COC(=O)C1C2CCC(OC1OC(=O)C1=CC=CC=C1)N2C', 0.8225806451612903), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8225806451612903), ('COC(=O)C1C2C(C)CC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8253968253968254), ('COC(=O)C12C(OC(=O)C3=CC=CC=C3)CC13CCC2N3C', 0.8253968253968254), ('COC(=O)C(C(C)OC(=O)C1=CC=CC=C1)C1CCCN1C', 0.8253968253968254), ('COC(=O)C1C2CCC(=CC1OC(=O)C1=CC=CC=C1)N2C', 0.8253968253968254), ('COC(=O)C1C2CCC(CC1OCC1=CC=CC=C1)N2C', 0.8253968253968254), ('CNC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8253968253968254), ('COC(=O)C1C2CCC(CC1NOC(=O)C1=CC=CC=C1)N2C', 0.828125), ('COC(=O)C1C2CCC(CC1(O)OC(=O)C1=CC=CC=C1)N2C', 0.828125), ('C=C(OC)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.828125), ('COC(=O)C1NC2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.8307692307692307), ('COC(=O)C1C2CCC(OCC1OC(=O)C1=CC=CC=C1)N2C', 0.8307692307692307), ('COC(=O)C1C(OC(=O)C2=CC=CC=C2)CC2CCC1(C)N2C', 0.8307692307692307), ('COC12CCC(CC(OC(=O)C3=CC=CC=C3)C1C=O)N2C', 0.8333333333333334), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CNC=C1)N2C', 0.8382352941176471), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CCC=C1)N2C', 0.8484848484848485), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC3=CC=C31)N2C', 0.8507462686567164), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC(O)CC=C1)N2C', 0.8529411764705882), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CN=CC=C1)N2C', 0.855072463768116), ('COC(=O)C1C2CC3C4=CC=C(C=C4)C(=O)OC1CC3N2C', 0.859375), ('COC(=O)C1C2CCC3(CC1OC(=O)C1=CC=C3C=C1)N2C', 0.8636363636363636), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC1)N2C', 0.8656716417910448), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1C)N2C', 0.8656716417910448), ('COC(=O)C1C2CCN(C)C1CCC1=CC=CC=C1C(=O)O2', 0.8656716417910448), ('COC(=O)C1C2CCC(CC1OC(=O)C1=C=CC=C1)N2C', 0.8676470588235294), ('COC(=O)C1C2CCC(CC1OC(=O)C1=NC=CC=C1)N2C', 0.8695652173913043), ('COC(=O)C1C2=CC=CC(=C2)C(=O)OCCC2CCC1N2C', 0.8695652173913043), ('COC(=O)C1C2CCC3CC1OC(=O)C1=C(C=CC=C1)CN32', 0.8695652173913043), ('COC(=O)C1C2CCC(CC13OC(=O)C1=CC=CC=C13)N2C', 0.8695652173913043), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=N1)N2C', 0.8695652173913043), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC(C)C1)N2C', 0.8714285714285714), ('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CCC1C)N2C', 0.8714285714285714), ('COC(=O)C(C1CCCN1C)C1CC2=CC=CC(=C2)C(=O)O1', 0.875), ('COC(=O)C1C2CCC(CC1OC1=CC=C(C=O)C=C1)N2C', 0.8787878787878788), ('COC(=O)C1C2CCC(CC1OC(O)C1=CC=CC=C1)N2C', 0.8787878787878788), ('CN1C2CCC1C1C(=O)OCC3=CC=CC(=C3)C(=O)OC1C2', 0.8840579710144928), ('COC(=O)C1C2CCC(CC1OC(=O)C1C=CC=C1)N2C', 0.90625), ('CC=CC=C(C)C(=O)OC1CC2CCC(C1C(=O)OC)N2C', 0.9285714285714286), ('C=CC=CC=CC(=O)OC1CC2CCC(C1C(=O)OC)N2C', 0.9285714285714286)]

When the list of candidates is populated and sorted, we need to choose the morphs that
will form the next :term:`generation <morph generation>` (the next leaves of the tree from which new morphs can be generated).
Here is an example implementation of a very simple filtering procedure:

..  code-block:: python
    :caption: A simple morph filter that selects only the first three closest morphs from the list.
    :name: filtering-morphs
    :linenos:

    # print the current candidates mask (all positions are on by default)
    print(tree.candidates_mask)

    # accept only the first three morphs in the sorted list (those with the lowest distance to target)
    mask = [False for x in tree.candidates_mask]
    mask[0] = True
    mask[1] = True
    mask[2] = True

    # save the new mask to the tree
    tree.candidates_mask = mask

    # show results
    print(tree.candidates_mask)
    print(
        [
            (x.smiles, x.dist_to_target)
            for idx,x in enumerate(tree.candidates)
            if tree.candidates_mask[idx] # get accepted molecules only
        ]
    )

Output:

..  code-block:: none

    (True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True)
    (True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False)
    [('COC(=O)C1C2CCC(CC1OCC1=CC=CC=C1)N2C', 0.5), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.7704918032786885)]

In :numref:`filtering-morphs` `candidates_mask` member can be easily changed by writing
a `list` or a `tuple` of new values into it. Here we simply select the first three morphs as the new :term:`morph generation`.
Note that the list of selected morphs is sorted from the molecule closest to the target to the farthest,
because we called the `sortMorphs()` method previously.

..  note:: The new mask must be the same length as the `candidates` member. If this requirement
        is not satisfied, an instance of `RuntimeError` is raised.

..  warning:: The mask should only be set after the morphs are sorted. If the mask is set and
        the order of morphs is changed, the mask will stay the same and will have to be updated
        to follow the new order.

..  seealso:: The library implements a few built-in filters. You can use the
        `filterMorphs()` method to invoke them. See the method's documentation for more information
        on the available filtering options.

.. _extend-prune:

Extending and Pruning
~~~~~~~~~~~~~~~~~~~~~

When we have the morphs we want to attach to the tree selected, we can call `extend()`
which will connect them to their respective parents making them the new leaves
of our tree:

..  code-block:: python
    :caption: Extending and pruning the tree.
    :name: extending-tree
    :linenos:

    # get the number of generations before
    print(tree.generation_count)

    tree.extend() # connect the accepted morphs to the tree as new leaves
    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda item : item[1]
        )
    )

    # get the number of generations after
    print(tree.generation_count)

    # check if a path was found
    print(tree.path_found)

    # run the pruning operation on the updated tree
    tree.prune()

Output:

..  code-block:: none

    [('COC(=O)C1C2CCC(CC1OCC1=CC=CC=C1)N2C', 0.5), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('CCOC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C', 0.7704918032786885)]
    1
    False

We can see that after extending the tree, the selected morphs (see :numref:`filtering-morphs`)
had become the new leaves and that the
:term:`morph generation` counter (`generation_count`) was set to one.

Because a tree generated in this way grows exponentially, a pruning strategy is needed in order
to keep the number of explored putative paths to a minimum by discarding those that are not getting any
closer to the :term:`target molecule`.

We call the molecule that have not generated any morphs closer to the target than itself
a :term:`non-producing molecule` and we can set
the number of generations to wait before removing its descendents
with the `non_producing_survive` parameter.

Tree pruning can be requested anytime by calling the `prune()` method. In our example, the method didn't prune
any paths, because the `non_producing_survive` parameter is set to 2 generations in this particular instance.

..  seealso:: In addition to the `non_producing_survive` parameter, there is the `max_morphs_total` parameter,
        which imposes a restriction on the maximum number of
        descendents of one :term:`non-producing molecule`. If the number of descendents
        reaches this threshold, the molecule is removed along with the descendents.

..  hint:: We could now repeat the process of generating, filtering, extending and pruning
        and have Molpher iteratively explore further in :term:`chemical space`. We could also read the `path_found` member
        at every iteration to check for the presence of the :term:`target molecule` in the tree and terminate the process,
        if a it is present. This would give us one complete implementation of a chemical space exploration algorithm.

Now we know everything that is needed to implement a chemical space exploration algorithm with
the *Molpher-lib* library. In the following sections, we describe more advanced topics
and introduce some built-in features of the library that can help to make some more complex tasks
(such as tree serialization) easier.

..  _operations:

Tree Operations
---------------

We call every action that is performed on an :term:`exploration tree` a *tree operation*.
This concept is represented in the library with the :class:`~molpher.core.operations.TreeOperation` abstract class and it
becomes useful when we run into a situation where we need to build
several exploration trees at once, want to reuse some existing code or store some interim results
of an ongoing exploration.

We can run any operation on a tree simply by supplying it to the
`runOperation()` method of the :py:class:`~molpher.swig_wrappers.core.ExplorationTree` class. Here is how to implement
the same workflow as in the preceding sections using operations:

..  code-block:: python
    :caption: Using operations to implement simple chemical space exploration.
    :emphasize-lines: 3-23,31
    :name: operations-example
    :linenos:

    from molpher.core.operations import *

    class MyFilterMorphs(TreeOperation):
        """
        A custom tree operation that accepts
        only the first three morphs
        (those with the lowest distance to target).

        """

        def __call__(self):
            """
            This method is called automatically by the tree.
            The tree this operation is being run on is accessible
            from `self.tree`.

            """

            mask = [False for x in self.tree.candidates_mask]
            mask[0] = True
            mask[1] = True
            mask[2] = True
            self.tree.candidates_mask = mask

    tree = ETree.create(source=cocaine, target=procaine) # create the tree

    # this list of tree operations defines one iteration
    iteration = [
        GenerateMorphsOper()
        , SortMorphsOper()
        , MyFilterMorphs()
        , ExtendTreeOper()
        , PruneTreeOper()
    ]

    # apply the operations in the list one by one
    for oper in iteration:
        tree.runOperation(oper)

    # observe the results
    print(tree.generation_count)
    print(tree.path_found)
    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
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

..  note:: Because the morphing algorithm is not deterministic and we initilized a new tree,
        the set of obtained morphs is slightly different from the one in the previous examples.

Most of the operations in :numref:`operations-example` are built-in operations (discussed below), but
we chose to define our own operation for the filtering step
(see the highlighted lines). We simply created a subclass of the :class:`~molpher.core.operations.TreeOperation`
abstract class and we overrode its
:py:meth:`~molpher.swig_wrappers.core.TreeOperation.__call__` method with the implementation we want.

Each operation can have a tree associated with it, but it is not necessary.
We can verify if a tree is associated with an operation by calling
its :meth:`~operations.TreeOperation.TreeOperation.getTree()`
method or accessing the `TreeOperation.tree` attribute of the class.
If there is no tree associated with the instance, they both return `None`.

..  note:: The built-in operations will raise a `RuntimeError`, if invoked without a tree attached to them.

Built-in Operations
~~~~~~~~~~~~~~~~~~~

A few operations are already defined in the library:

    - `GenerateMorphsOper`
    - :py:class:`~operations.SortMorphsOper.SortMorphsOper`
    - `FilterMorphsOper`
    - `FindLeavesOper`
    - `ExtendTreeOper`
    - `PruneTreeOper`
    - `TraverseOper`

They are all dervied from :class:`~molpher.swig_wrappers.core.TreeOperation` and contain
the full set of operations performed on a tree in
the original Molpher algorithm as published in [1]_. Therefore, the original algorithm can be
implemented using those operations.

In the next part of the tutorial, we will pay particular attention to the `TraverseOper` operation. It differs
from the others, because it uses a callback function to perform actions on molecules
in the tree and is, therefore, very useful for debugging and saving exporting various data (see `tree-traversal`).

For more details on the other operations, see the designated pages in the documentation.

.. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

..  _tree-traversal:

Traversing the Tree
^^^^^^^^^^^^^^^^^^^

A special place among the operations belongs to the `TraverseOper` class. It does not directly implement a part
of a morphing algorithm, but serves as a means of traversing molecules in a tree and reading/modifying them
as needed:

..  code-block:: python
    :caption: Traversing the tree using a callback.
    :name: traverse-example
    :linenos:

    class MyCallback(TraverseCallback):
        """
        This callback just prints some information
        about the molecules in the tree.

        """

        def __call__(self, morph):
            """
            Method called on each morph in the tree
            -- starting from the root to leaves.

            """

            if not morph.getParentSMILES():
                print("# Root #")
            else:
                print('# Morph #')
                print('Parent:', morph.getParentSMILES())
            print('SMILES: ', morph.getSMILES())
            print('Descendents: ', morph.getDescendants())

    callback = MyCallback() # initialize a callback
    traverse = TraverseOper(callback) # attach it to a tree traversal operation
    tree.runOperation(traverse) # run the operation

Output:

..  code-block:: none

    # Root #
    SMILES:  COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ('COC(=O)C1C(OC(=O)C2=CC=CC=C2)CC2CCN1N2C', 'COCC1C2CC(=O)C(CC1OC(=O)C1=CC=CC=C1)N2C', 'COCC1C2CCC(C(=O)C1OC(=O)C1=CC=CC=C1)N2C')
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COCC1C2CCC(C(=O)C1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COCC1C2CC(=O)C(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)C1C(OC(=O)C2=CC=CC=C2)CC2CCN1N2C
    Descendents:  ()

In :numref:`traverse-example` we derive from the `TraverseCallback` class, an abstract class with
an abstract method :meth:`~molpher.swig_wrappers.core.TraverseCallback.__call__`. This method takes one argument,
which is a :class:`~molpher.core.MolpherMol.MolpherMol` instance
of a molecule in the tree. We need to override this method in our derived class in order to implement our own
behaviour.

The callback is then associated with a `TraverseOper` instance, which can be run on a tree as any other
tree operation. When the operation is run it traverses the tree from the root to the leaves and injects
every molecule it encounters into our implementation of the
:meth:`~molpher.swig_wrappers.core.TraverseCallback.__call__` method.

..  note:: We can also pass a SMILES string to the `TraverseOper` constructor. In that case, a subtree will be traversed
        using the specified molecule as the root of the subtree.

There is also a much more convenient way to traverse the tree. Because, the `ExplorationTree` class implements
the `traverse()` method, we can simply take any python callable and use it instead of the
:meth:`~molpher.swig_wrappers.core.TraverseCallback.__call__` method.
However, under the hood it does the same thing as we did in :numref:`traverse-example`.
Therefore, the above code can be turned into:

..  code-block:: python
    :caption: Traversing the tree using a callback -- the simple version.
    :name: short-traverse-example
    :linenos:

    def process(morph):
        """
        Prints some information
        about the molecules in the tree.

        """

        if not morph.getParentSMILES():
            print("# Root #")
        else:
            print('# Morph #')
            print('Parent:', morph.getParentSMILES())
        print('SMILES: ', morph.getSMILES())
        print('Descendents: ', morph.getDescendants())

    tree.traverse(process) # use the traverse method to run the callback function

Output:

..  code-block:: none

    # Root #
    SMILES:  COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ('COC(=O)C1C(OC(=O)C2=CC=CC=C2)CC2CCN1N2C', 'COCC1C2CC(=O)C(CC1OC(=O)C1=CC=CC=C1)N2C', 'COCC1C2CCC(C(=O)C1OC(=O)C1=CC=CC=C1)N2C')
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COCC1C2CCC(C(=O)C1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COCC1C2CC(=O)C(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)C1C(OC(=O)C2=CC=CC=C2)CC2CCN1N2C
    Descendents:  ()

..  _templates-snapshots:

Tree Templates and Snapshots
----------------------------

We don't always have to initialize :term:`morphing parameters` by hand. We can use an :term:`XML template` instead.
Here is an example of a template file (you can download this one
from :download:`here <../../../src/python/molpher/examples/cocaine-procaine-template.xml>`):

..  literalinclude:: ../../../src/python/molpher/examples/cocaine-procaine-template.xml
    :language: xml
    :caption: A complete XML template file.
    :name: template-file
    :linenos:

An :term:`XML template` is similar to a configuration file and can be loaded
just like a snapshot (see :numref:`loading-snapshot`), but the resulting tree
will only contain the :term:`source molecule` as its root.

..  code-block:: python
    :caption: Loading a template and saving a built tree as a XML snapshot.
    :name: saving-snapshot
    :linenos:

    template_file = 'cocaine-procaine-template.xml'

    tree = ETree.create(template_file)
    print(tree.params)

    # apply the tree operations
    for oper in iteration:
        tree.runOperation(oper)

    print(
        sorted( # grab the new leaves as a list sorted according to their distance from target
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda x : x[1]
        )
    )

    # save the tree in a snapshot file
    tree.save('snapshot.xml')

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
        'source': 'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C',
        'target': 'CCN(CC)CCOC(=O)C1=CC=C(N)C=C1',
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

In the above example we loaded an :term:`XML template`, created a tree from it, extended the tree and
serialized it as a snapshot. We can see that all the parameters are the same as in the :term:`XML template` and that
the resulting tree can be built using the same list of operations
as in :numref:`operations-example`.

In :numref:`saving-snapshot` we also serialized our tree instance to disk
with the :py:meth:`~molpher.swig_wrappers.core.ExplorationTree.save` method.
The saved tree can be later reconstructed with the
:meth:`~molpher.core.ExplorationTree.ExplorationTree.create` factory method:

..  code-block:: python
    :caption: Loading a snapshot of a previously generated tree.
    :name: loading-snapshot
    :linenos:

    new_tree = ETree.create('snapshot.xml') # create a new tree from the saved snapshot
    print(new_tree.params)
    print(
        sorted( # grab the leaves in the created tree (these should be the same as those in the original tree)
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
        'source': 'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C',
        'target': 'CCN(CC)CCOC(=O)C1=CC=C(N)C=C1',
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

Now we wrap up this tutorial with two simple chemical space exploration implementations
(see :numref:`complete-example`
and :numref:`bidirectional-example`).

..  _simple:

Using the Tutorial to Implement a Search Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the tutorial we showed how to implement
one step of a simple algorithm that searches
for a path in :term:`chemical space` between *cocaine* and *procaine*.
Transforming the code into a full exploration algorithm is pretty much straightforward:

..  literalinclude:: ../../../src/python/molpher/examples/experiment.py
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

The above implementation is nothing more than just the tutorial code bits inside a loop.
The loop checks if a path was found at each iteration.
If the path is found, it backtracks through the tree
and prints out a sequence of molecules lying on the path.

..  _bidirectional:

Implementing a Bidirectional Search Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second example we have here is a little bit more elaborate,
but implements a very simple idea. Instead of one exploration tree,
we build two trees that each search for a path to the closest molecule in the other:

..  literalinclude:: ../../../src/python/molpher/examples/bidirectional.py
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

This bidirectional algorithm uses the built-in operations to facilitate the search,
but does one extra procedure after
an iteration is completed -- it changes the target molecules of the trees.
When the new leaves are connected both trees are traversed and molecules
closest to the current target are identified in each. The closest molecule from one tree is then
set as the new target for the tree searching in the opposite direction and vice versa.

In :numref:`bidirectional-example` we also use the `time.clock` function to measure the execution
times of each potentially time-consuming operation.

Summary
-------

If you have read the tutorial all the way to here, you now probably have a decent idea on what the library does
and how to use it. If you have any suggestions on how to improve it or bug reports, please send them to
sichom@vscht.cz or create an issue on the issue tracker. All help on the project is much appreciated.

..  todo:: link to issue tracker

