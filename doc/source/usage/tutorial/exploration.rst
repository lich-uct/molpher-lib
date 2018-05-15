..  _exploration:

Implementing Chemical Space Exploration with Molecular Morphing
===============================================================

..  todo:: divide this doc into multiple smaller ones and use TOC to show contents

.. contents:: Table of Contents

In the :doc:`previous tutorial <morphing>`, we introduced some of the low level
Molpher-lib features that enable its user to modify molecular
structures in a randomized manner. In this part of the usage tutorial, we focus on
high level features that handle chemical space exploration on
larger scale by iterative application of :term:`chemical operators`.

The core data structure responsible for building and maintaining suggested
:term:`chemical space paths <chemical space path>` is the :term:`exploration tree`.
The exploration tree contains all generated morphs and is rooted at the :term:`source molecule`
and is grown by iteratively applying chemical operators
first on the source molecule, then on its descendants and so on. This way
multiple chemical space paths can be created and evaluated.

As the paths are generated, the compounds on them can be characterized by
the value of an objective function. This function should somehow formalize
the relationship between a chemical structure
and its fitness for the task at hand. Therefore,
interesting paths can be prioritized over others.
This way, the exploration tree can hopefully lead us towards interesting new compounds
without having to enumerate all possible chemistry.

In the original Molpher algorithm, the objective function is
simply the structural distance between morphs and a :term:`target molecule`.
This way Molpher is able to connect pairs of molecules by iteratively generating
morphs from the source and prioritizing those that are getting closer to the target.
As a result, the algorithm converges to the target structure and stops
once the target is identified among the generated morphs (:numref:`morphing-explore`). This
algorithm is called the *classic* algorithm and is
implemented in the `molpher.algorithms.classic` module.

..  seealso:: `alg-classic`

..      figure:: morphing.png
        :scale: 100%
        :name: morphing-explore

        Schematic depiction of the original algorithm published by Hoksza et al. [1]_.
        New :term:`morphs <morph>` are generated with `chemical operators`
        applied to the :term:`source molecule` (M\ :sub:`S`\ ) until the :term:`target
        molecule` (M\ :sub:`T`\) is found.

..  seealso:: Examples in this tutorial are included in a Jupyter Notebook
    which can be downloaded :download:`here <../../../notebooks/exploration.ipynb>`
    viewed `here <../../_static/exploration.html>`_. The example SDF can be found
    :download:`here <../../../notebooks/captopril.sdf>`.

.. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

..  _tree-create:

Creating an Exploration Tree and Setting Morphing Parameters
------------------------------------------------------------

.. py:currentmodule:: molpher.core

In Molpher-lib, the tree is implemented as the `molpher.core.ExplorationTree` class.
Following the captopril example from the previous tutorial,
the simplest way to generate an exploration tree instance is
to call the :meth:`~.core.ExplorationTree.create` factory method:

..  code-block:: python
    :linenos:

    from molpher.core import ExplorationTree as ETree
    from molpher.core import MolpherMol

    captopril = MolpherMol("captopril.sdf")
    tree = ETree.create(source=captopril)

This code simply initializes the tree from the supplied :class:`MolpherMol` instance.
At the moment the tree is pretty simple. It only contains the captopril structure as its
current leaves:

..  code-block:: python

    tree.leaves[0].asRDMol()

Output:

..  figure:: captopril_numbered.png

We can manipulate this tree and read data from it. Let's start by printing out the
:term:`source molecule`:

..  code-block:: python

    print('Source: ', tree.params['source'])

Output:

..  code-block:: none

    Source:  CC(CS)C(=O)N1CCCC1C(=O)O

The `params` member is a dictionary used to set and store `morphing parameters`:

..  code-block:: python

    print(tree.params)

Output:

..  code-block:: python

    {
        'source': 'CC(CS)C(=O)N1CCCC1C(=O)O',
        'target': None,
        'operators': (
            'OP_ADD_ATOM',
            'OP_REMOVE_ATOM',
            'OP_ADD_BOND',
            'OP_REMOVE_BOND',
            'OP_MUTATE_ATOM',
            'OP_INTERLAY_ATOM',
            'OP_BOND_REROUTE',
            'OP_BOND_CONTRACTION'
        ),
        'fingerprint': 'FP_MORGAN',
        'similarity': 'SC_TANIMOTO',
        'weight_min': 0.0,
        'weight_max': 100000.0,
        'accept_min': 50,
        'accept_max': 100,
        'far_produce': 80,
        'close_produce': 150,
        'far_close_threshold': 0.15,
        'max_morphs_total': 1500,
        'non_producing_survive': 5
    }

As we can see there is quite a lot of parameters that we can set,
but most of these affect the exploration process only if
some parts of the library are used in the context of the tree, especially tree operations
which we will discuss :ref:`later <operations>`. The most important parameters
will be explained in this tutorial, but you can see the
documentation for the :py:class:`~ExplorationData.ExplorationData` class
(especially :numref:`param-table`) for a more detailed reference.

We can adjust the morphing parameters during runtime as we like.
All we need to overwrite the `params` attribute
of our tree instance with a new dictionary:

..  code-block:: python

    # change selected parameters using a dictionary
    tree.params = {
        'non_producing_survive' : 2
        , 'weight_max' : 500.0
    }
    print(tree.params)

Output:

..  code-block:: python

    {
        'source': 'CC(CS)C(=O)N1CCCC1C(=O)O',
        'target': None,
        'operators': (
            'OP_ADD_ATOM',
            'OP_REMOVE_ATOM',
            'OP_ADD_BOND',
            'OP_REMOVE_BOND',
            'OP_MUTATE_ATOM',
            'OP_INTERLAY_ATOM',
            'OP_BOND_REROUTE',
            'OP_BOND_CONTRACTION'
        ),
        'fingerprint': 'FP_MORGAN',
        'similarity': 'SC_TANIMOTO',
        'weight_min': 0.0,
        'weight_max': 500.0,
        'accept_min': 50,
        'accept_max': 100,
        'far_produce': 80,
        'close_produce': 150,
        'far_close_threshold': 0.15,
        'max_morphs_total': 1500,
        'non_producing_survive': 2
    }

Here we just tightened the constraints on molecular weight
for the morphs that we allow to be incorporated in the tree
(applied only if the :class:`~.operations.FilterMorphsOper.FilterMorphsOper` operation or `filterMorphs` method
are used with certain options) and we decreased the number of acceptable 'non-producing'
:term:`morph generations <morph generation>` to 2. Non-producing generations are
generations of morphs that has not improved in the objective function
(e.g. structural distance). See :numref:`param-table` for details.
One thing to note is that if we supply an incomplete set of parameters
(like in the example above), only the parameters specified in the
supplied dictionary will be changed. Parameters not mentioned in
this dictionary will remain the same as before the assignment.

..  warning:: Changing individual values in the `params` dictionary will have no effect.
    You always need to store a dictionary instance in it. This is because the value
    is regenerated every time the attribute is accessed to always reflect the current
    set of parameters valid for the current instance.

..  seealso:: :py:class:`~ExplorationData.ExplorationData`

Generating Morphs and Extending the Exploration Tree
----------------------------------------------------

This part of the tutorial outlines the steps
involved in one iteration of a possible exploration algorithm.
We explain how to generate new morphs from the leaves of the tree,
how they can be filtered and how the tree is
extended by attaching the chosen morphs as the
next generation of leaves. We also show how the unfavorable paths (or their parts)
can later be removed from the growing tree.

Generating and Manipulating Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Previously <tree-create>`, we showed how to initialize an exploration tree.
Now that we have one, we can take a look at how to use it for :term:`chemical space`
exploration.

Let us generate a few :term:`morphs <morph>` from the current leaves of the tree.
Currently, the tree has just one leaf (our source molecule, captopril):

..  code-block:: python

    print(tree.leaves) # show the current leaves of the tree (only the source so far)
    tree.leaves[0].asRDMol()

Output:

..  code-block:: none

    (<molpher.core.MolpherMol.MolpherMol at 0x7f37a8a6d930>,)

..  figure:: captopril.png

Since we already have the built-in operators available by default,
we can generate new structures from this starting molecule like so:

..  code-block:: python
    :caption: Generating candidate morphs with an exploration tree.
    :name: gen-morphs-simple
    :linenos:

    tree.generateMorphs()
    print(len(tree.candidates))

Output:

..  code-block:: none

    28

The :meth:`~.core.ExplorationTree.ExplorationTree.generateMorphs()` method tells
the tree to generate some :term:`morphs <morph>`
from the current leaves for us. The number of generated morphs
will depend mostly on the `far_produce`, `close_produce` and `far_close_threshold`
parameters. However, it also depends on other factors. For example, some structures
might not be parsed correctly and, thus, might not make it to the final list.
Also, a different number of morphs can be generated each time the method is run. That si due to
the non-deterministic character of the morphing algorithm which chooses the morphing operators to
use and parts of the structure to modify randomly. Duplicate
molecules (based on the canonical smiles string) are also removed.

We can access the newly generated morphs from the `candidates`
member of the tree instance. It is a `tuple` of :py:class:`~MolpherMol.MolpherMol` instances.
These instances can be used to read and manipulate the generated morphs or
the compounds currently present in the tree.

..  attention:: The molecules saved in the `candidates` attribute of the tree actually do not
    belong to the tree just yet. See :ref:`extend-prune` for more information on
    how tree ownership is assigned to molecules.

Sorting and Filtering Morphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The order of the newly generated molecules in the `candidates` list has a meaning
for the search algorithm. The higher the position of a morph in this list,
the bigger the probability that we will attach it to the tree as a new leaf
(see `PROBABILITY` for details). Therefore, by sorting this list according
to a given objective function, it is possible to push the algorithm
into convergence.

As of yet, the only way to sort the
generated morphs is by calling the :py:meth:`~molpher.swig_wrappers.core.ExplorationTree.sortMorphs`
method on the tree instance or using the :py:class:`~operations.SortMorphsOper.SortMorphsOper` operation
(see :ref:`operations` for more). This sorts the molecules in the order of increasing value
of the objective function. Right now, the objective function value used
by the `sortMorphs` operation is the value of the `dist_to_target` property
of the given :class:`~.core.MolpherMol` instance. By default, the value
of this property indicates the structural distance of the morph in question
from the target molecule.

Since in our case we did not specify a target, the `dist_to_target` property will
be set to the same value for all molecules:

..  code-block:: python

    {x.dist_to_target for x in tree.candidates}

Output:

..  code-block:: none

    {1.7976931348623157e+308}

However, the `dist_to_target` property value can be changed to basically any float number.
If users have a custom objective function they want to use in the search,
it is possible to write their calculated value to this property and as long as lower value means
better fitness of the morph, it should still make sense to use this function
in the context of Molpher-lib and its data structures.

..  note:: This part of the interface (especially the naming convention) will
    likely change in the future since a more general approach is needed in order to make sorting more
    customizable.

:ref:`Previously <morphing-bulk>`, we introduced the concept of morph collectors,
which are special functions that can be used to intercept morphs as they are created.
Setting the objective function value for each morph is a perfect use case for them:

..  code-block:: python
    :linenos:

    def sascore_as_obj(morph, operator):
        morph.dist_to_target = morph.sascore

    tree.generateMorphs([sascore_as_obj])
    print(len(tree.candidates))
    {x.dist_to_target for x in tree.candidates}

Output (some results omitted):

..  code-block:: none

    28

    [4.304767951403637,
     4.161464339486345,
     3.8871106534610247,
     ...
     4.336516757110866,
     4.187278867161213,
     3.8893996483733346]

This code is essentially the same as in :numref:`gen-morphs-simple`,
but this time we chose to add a custom morph collector
(the :code:`sascore_as_obj` function),
which will set a molecule's :term:`SAScore` as the objective
function value for all generated morphs. Notice
that the :meth:`~.core.ExplorationTree.generateMorphs` method takes a list of collectors
so it is possible to chain them together.
They are applied in the order of appearance in the list.

As you can see, the list of candidates is not sorted, yet.
We need to call the :meth:`~.core.ExplorationTree.sortMorphs` method to do that:

..  code-block:: python
    :caption: Sorting morphs according to the value of the objective function.
    :linenos:

    tree.sortMorphs()

    [
        (x.smiles, x.dist_to_target)
        for idx,x in enumerate(tree.candidates)
    ]

Output (some results omitted):

..  code-block:: none

    [('CC(C)C(=O)N1CCCC1C(=O)O', 3.3191105796423788),
     ('CC(CO)C(=O)N1CCCC1C(=O)O', 3.6043596148886445),
     ('CNC(CS)C(=O)N1CCCC1C(=O)O', 3.7075484704945465),
    ...
     ('O=C(O)C1CCCN1C(=O)C(CS)CBr', 4.404706288979395),
     ('CC(CSI)C(=O)N1CCCC1C(=O)O', 4.412717102918789),
     ('O=C(O)C1CCCN1C(=O)C(CF)CS', 4.420781866555153)]

Therefore, now the list of candidates is sorted according
to their synthetic accessibility (compounds
that are easier to prepare *in vitro* should have lower scores).

Now, we need to choose the morphs that
will form the next :term:`generation <morph generation>`.
The `candidates_mask` property of :class:`~.core.ExplorationTree`
serves exactly this purpose. Each position in this list corresponds to
one molecule in `candidates` and indicates
whether this molecule should be considered when
attaching new leaves to the tree (`True`) or not (`False`).
Here is an example implementation of a very simple filtering procedure:

..  code-block:: python
    :caption: A simple morph filter that selects only the first three closest morphs from the list.
    :name: filtering-morphs
    :linenos:

    # print the current candidates mask (all positions are on by default)
    print("Old mask:", tree.candidates_mask)

    # accept only the first ten morphs in the sorted list (those with the lowest distance to target)
    new_mask = [True if idx < 10 else False for idx, x in enumerate(tree.candidates_mask)]

    # save the new mask to the tree
    tree.candidates_mask = new_mask

    # show results
    print("New mask:", tree.candidates_mask)
    print("Molecules that passed the filter:")
    [
        (x.smiles, x.dist_to_target)
        for idx,x in enumerate(tree.candidates)
        if tree.candidates_mask[idx] # get molecules that passed the filter only
    ]

Output:

..  code-block:: none

    Old mask: (True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True)
    New mask: (True, True, True, True, True, True, True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False)
    Molecules that passed the filter:

    [('CC(C)C(=O)N1CCCC1C(=O)O', 3.3191105796423788),
     ('CCC(C)C(=O)N1CCCC1C(=O)O', 3.404002369297247),
     ('CSCCC(=O)N1CCCC1C(=O)O', 3.613205289055311),
     ('CC(CS)C(=O)N1CCCC1C(=O)O', 3.804751376555311),
     ('O=C(O)C1CCCN1C(=O)CCS', 3.8871106534610247),
     ('O=C(O)C1CCCN1C(=O)CCCS', 3.9220880467166013),
     ('CC(S)CC(=O)N1CCCC1C(=O)O', 3.9366697951036973),
     ('O=C(O)C1CCCN1C(=O)C1CSC1', 3.9784865729838823),
     ('CC(S)C(=O)N1CCCC1C(=O)O', 3.9938638051851627),
     ('CC(NCS)C(=O)N1CCCC1C(=O)O', 4.076862613435724)]

In :numref:`filtering-morphs`, `candidates_mask` member was changed by writing
a `list` or a `tuple` of new values into it. This way we were able to select only the first ten morphs
that have the best `SAScore` value.

..  warning:: The mask should only be set after the morphs are sorted. If the mask is set and
    the order of morphs is changed, the mask will stay the same and will have to be updated
    to follow the new order.

..  seealso:: The library implements a few built-in filters. You can use the
    `filterMorphs()` method or :class:`~.operations.FilterMorphsOper` operation to invoke them.
    See the method's documentation for more information
    on the available filtering options.

.. _extend:

Extending
~~~~~~~~~

When we have the morphs selected, we can call
the `extend()` method. This will connect them to their respective parents
in our tree and they will become a new set of leaves:

..  code-block:: python
    :caption: Extending the exploration tree with new morphs.
    :name: extending-tree
    :linenos:

    # get the number of generations before
    print(tree.generation_count)

    tree.extend() # connect the accepted morphs to the tree as new leaves

    # get the number of generations after
    print(tree.generation_count)

    # grab the new leaves as a list sorted according to their distance from target
    sorted(
        [
            (x.getSMILES(), x.getDistToTarget())
            for x in tree.leaves
        ], key=lambda item : item[1]
    )

Output:

..  code-block:: none

    WARNING: Candidate morph: CC(CS)C(=O)N1CCCC1C(=O)O already present in the tree. Skipping...
    0
    1

    [('CC(C)C(=O)N1CCCC1C(=O)O', 3.3191105796423788),
     ('CCC(C)C(=O)N1CCCC1C(=O)O', 3.404002369297247),
     ('CSCCC(=O)N1CCCC1C(=O)O', 3.613205289055311),
     ('CSC(C)C(=O)N1CCCC1C(=O)O', 3.8501001628456333),
     ('O=C(O)C1CCCN1C(=O)CCS', 3.8871106534610247),
     ('CCC(CS)C(=O)N1CCCC1C(=O)O', 3.8893996483733346),
     ('CSCC(C)C(=O)N1CCCC1C(=O)O', 3.916140148729842),
     ('O=C(O)C1CCCN1C(=O)CCCS', 3.9220880467166013),
     ('CC(S)CC(=O)N1CCCC1C(=O)O', 3.9366697951036973)]

We can see that after extending the tree, the selected morphs (see :numref:`filtering-morphs`)
had become the new leaves and that the tree's
:term:`morph generation` counter (`generation_count`) was increased by one. We also
got a warning about one structure not being added to the tree. It is the structure
of captopril itself, which is already there. Thus, it is automatically skipped to prevent us from going
in circles.

If we want to, we can generate an image depicting the new leaves and the operators used to create them like so:

..  code-block:: python

    from rdkit.Chem.Draw import MolsToGridImage

    def get_locked_atoms(mol):
        return [(idx, atm) for idx, atm in enumerate(mol.atoms) if atm.is_locked]

    def show_mol_grid(mols):
        locked_atoms = [[y[0] for y in get_locked_atoms(x)] for x in mols]
        return MolsToGridImage(
            [x.asRDMol() for x in mols]
            , subImgSize=(250,200)
            , highlightAtomLists=locked_atoms
            , legends=[x.parent_operator for x in mols]
        )

    show_mol_grid(tree.leaves)

Output:

..  figure:: leaves.png

Note that the generated morphs satisfy the locks placed on the signature -pril substructure
in the original SDF file. Therefore, the tree is guaranteed to only contain structures
that have this structural pattern.

Therefore, by iterative application of the commands above, we would be able to generate
many possible structures of novel -pril compounds. This could prove useful
while exploring the structure-activity relationship in the development
of new ACE inhibitors, for example.

.. _prune:

Pruning
~~~~~~~

However, there is also the problem that the tree will grow exponentially
if we keep adding new morphs in this way. Thus, we will need a strategy
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

    Generated 67 morphs.
    1
    False
    [('COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C', 0.7068965517241379), ('COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C', 0.7142857142857143), ('COC(=O)C1C(COC(=O)C2=CC=CC=C2)C2CCC1N2C', 0.7586206896551724)]

..  note:: Because the morphing algorithm is not deterministic and we initilized a new tree,
    the set of obtained morphs is different from the one in the previous examples.

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

    - :py:class:`~operations.GenerateMorphsOper.GenerateMorphsOper`
    - :py:class:`~operations.SortMorphsOper.SortMorphsOper`
    - :py:class:`~operations.FilterMorphsOper.FilterMorphsOper`
    - :py:class:`~operations.FindLeavesOper.FindLeavesOper`
    - :py:class:`~operations.ExtendTreeOper.ExtendTreeOper`
    - :py:class:`~operations.PruneTreeOper.PruneTreeOper`
    - :py:class:`~operations.TraverseOper.TraverseOper`
    - :py:class:`~operations.CleanMorphsOper.CleanMorphsOper`

They are all dervied from :class:`~molpher.swig_wrappers.core.TreeOperation` and contain
the full set of operations performed on a tree in
the original Molpher algorithm as published in [1]_. Therefore, the original algorithm can be
implemented using those operations.

In the next part of the tutorial, we will pay particular attention to the
:py:class:`~operations.TraverseOper.TraverseOper` operation. It differs
from the others, because it uses a callback function to perform actions on molecules
in the tree and is, therefore, very useful for debugging and saving exporting various data (see `tree-traversal`).

For more details on the other operations, see the designated pages in the documentation.

..  _tree-traversal:

Traversing the Tree
^^^^^^^^^^^^^^^^^^^

A special place among the operations belongs to the :py:class:`~operations.TraverseOper.TraverseOper`
class. It does not directly implement a part
of a morphing algorithm, but serves as a means of traversing molecules in a tree and reading/modifying them
as needed:

..  code-block:: python
    :caption: Traversing the tree using a callback.
    :name: traverse-example
    :linenos:

    from molpher.core.operations.callbacks import TraverseCallback

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
    traverse = TraverseOper(callback=callback) # attach it to a tree traversal operation
    tree.runOperation(traverse) # run the operation

Output:

..  code-block:: none

    # Root #
    SMILES:  COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    Descendents:  ('COC(=O)C1C(COC(=O)C2=CC=CC=C2)C2CCC1N2C', 'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C', 'COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C')
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)C1C(COC(=O)C2=CC=CC=C2)C2CCC1N2C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C
    SMILES:  COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C
    Descendents:  ()
    Descendents:  ()

..  note:: The tree traversal algorithm uses multiple threads. Therefore,
    the output might not be synchronized.

In :numref:`traverse-example` we derive from the `TraverseCallback` class, an abstract class with
an abstract method :meth:`~molpher.swig_wrappers.core.TraverseCallback.__call__`. This method takes one argument,
which is a :class:`~molpher.core.MolpherMol.MolpherMol` instance
of a molecule in the tree. We need to override this method in our derived class in order to implement our own
behaviour.

The callback is then associated with a :py:class:`~operations.TraverseOper.TraverseOper`
instance, which can be run on a tree as any other
tree operation. When the operation is run it traverses the tree from the root to the leaves and injects
every molecule it encounters into our implementation of the
:meth:`~molpher.swig_wrappers.core.TraverseCallback.__call__` method.

..  note:: We can also pass a SMILES string to the :py:class:`~operations.TraverseOper.TraverseOper`
    constructor. In that case, a subtree will be traversed
    using the specified molecule as the root of the subtree.

There is also a much more convenient way to traverse the tree. Because, the `ExplorationTree` class implements
the :meth:`molpher.core.ExplorationTree.ExplorationTree.traverse()` method, we can simply take any python callable and use it instead of the
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
    Descendents:  ('COC(=O)C1C(COC(=O)C2=CC=CC=C2)C2CCC1N2C', 'COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C', 'COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C')
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)CC1CCC(CCOC(=O)C2=CC=CC=C2)N1C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)C1C2CCC(CC1OC(=O)C1=CC=C(N)C=C1)N2C
    Descendents:  ()
    # Morph #
    Parent: COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C
    SMILES:  COC(=O)C1C(COC(=O)C2=CC=CC=C2)C2CCC1N2C
    Descendents:  ()

..  _templates-snapshots:

Tree Templates and Snapshots
----------------------------

We don't always have to initialize :term:`morphing parameters` by hand. We can use an :term:`XML template` instead.
Here is an example of a template file (you can download this one
from :download:`here <../../../../src/python/molpher/examples/cocaine-procaine-template.xml>`):

..  literalinclude:: ../../../../src/python/molpher/examples/cocaine-procaine-template.xml
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

    The new iteration has been created from template:
    source: CN1[C@H]2CC[C@@H]1[C@@H](C(=O)OC)[C@@H](OC(=O)c1ccccc1)C2
    target: O=C(OCCN(CC)CC)c1ccc(N)cc1

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
    Generated 66 morphs.
    [('CN1C2CCC1C(C(=O)OCN)C(OC(=O)C1=CC=CC=C1)C2', 0.7777777777777778), ('CCN1C2CCC1C(C(=O)OC)C(OC(=O)C1=CC=CC=C1)C2', 0.7936507936507937), ('CN1C2CCC1C(C(=O)ON)C(OC(=O)C1=CC=CC=C1)C2', 0.8064516129032258)]

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
    [('CN1C2CCC1C(C(=O)OCN)C(OC(=O)C1=CC=CC=C1)C2', 0.7777777777777778), ('CCN1C2CCC1C(C(=O)OC)C(OC(=O)C1=CC=CC=C1)C2', 0.7936507936507937), ('CN1C2CCC1C(C(=O)ON)C(OC(=O)C1=CC=CC=C1)C2', 0.8064516129032258)]

Implementing a Chemical Space Exploration Algorithm
---------------------------------------------------

Previously, we showed how to implement
one step of a simple exploration algorithm.
Transforming the code into a full exploration algorithm is pretty much straightforward:

..  literalinclude:: ../../../../src/python/molpher/examples/algorithms_basics.py
        :language: python
        :caption: Example implementation of a pathfinding algorithm that searches
            for a path in :term:`chemical space` between *cocaine* and *procaine*.
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
If it was found, we backtrack through the tree
and print out a sequence of molecules on the path.

Summary
-------

If you read the tutorial all the way to here, you now probably have a decent idea on what the library does
and how to use it. If you have any suggestions on how to improve it or bug reports, please submit them to
the `issue tracker <https://github.com/lich-uct/molpher-lib/issues>`_. Any help on the project is much appreciated.
