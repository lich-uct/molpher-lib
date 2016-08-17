Example Exploration Algorithm Implementations
---------------------------------------------

The library also provides complete implementations of a few
exploration algorithms (located in :mod:`molpher.algorithms`). This section briefly describes these algorithms
and shows how to use them. It is also recommended to look at
the source code documentation of the :mod:`~molpher.algorithms` package,
which contains some generally useful modules
(such as :mod:`~molpher.algorithms.functions` and :mod:`~molpher.algorithms.operations`).

Classic (Original) Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the algorithm published with the :term:`Molpher`
program [1]_. It is the most basic algorithm
that can be implemented and a lot of fundamental concepts of
Molpher-lib are based on it. The other algorithms presented in this part of the documentation
are modifications or extensions of the original. It is, therefore,
recommended to read this section before moving on to the following ones.

The classic algorithm uses a single tree rooted at the `source molecule` (M\ :sub:`S`).
Using predefined `chemical operators`, it generates a new generation of structures from the leaves of this tree
at every iteration (see :numref:`morphing`).
The new list of candidates is then sorted (using the
:class:`~molpher.core.operations.SortMorphsOper.SortMorphsOper` operation)
in ascending order according to their structural distance
from the `target molecule` (M\ :sub:`T`).
Finally, the new structures are filtered with the
:class:`~molpher.core.operations.FilterMorphsOper.FilterMorphsOper`
operation and survivors are attached to their parents in the tree.
This is repeated until the target molecule is generated and appended to the tree.

As the computation proceeds, some tree branches are also pruned
(see :class:`~molpher.core.operations.PruneTreeOper.PruneTreeOper`). This is
necessary to prevent exponential growth of the tree by removing paths
that are not getting closer to the `target molecule`.

If you wish to know more about the steps involved
in this algorithm, a detailed description can be found in [1]_.

..      figure:: morphing.png
        :scale: 100%
        :name: morphing

        Schematic depiction of the original algorithm published by Hoksza et al. [1]_.
        New `morphs <morph>` are generated with the `chemical operators` until the target
        molecule is found.

The classic algorithm is available through the `molpher.algorithms.classic` package
and can be used as follows:

..  literalinclude:: ../../../src/python/molpher/examples/classic.py
        :language: python
        :caption: Usage example of the original algorithm.
        :name: classic-example
        :linenos:

The most essential ingredients in the example above are the :samp:`run()`
function and the appropriate :samp:`Settings` class.

Every algorithm in the library has a :samp:`run()` function
which calls the appropriate algorithm and saves
the output of the computation. The :samp:`run()` function
also supplies parameters to the tree and manages some 'global'
aspects when a path is generated (such as the maximum number of threads to use).

The behaviour of :samp:`run()` is configured using a :samp:`Settings`
class which is nothing more than a set of parameters wrapped into
a class object. Some algorithms might define their own :samp:`Settings`
classes with parameters specific to that algorithm. In this implementation, however,
the most general `Settings` class is used.

In this instance, :meth:`~molpher.algorithms.classic.run`
outputs a pickled `list` of SMILES strings
which represents the resulting path into `storage_dir`.


.. [1] Hoksza D., Škoda P., Voršilák M., Svozil D. (2014) Molpher: a software framework for systematic chemical space exploration. J Cheminform. 6:7.
        `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24655571>`_, `DOI <http://www.jcheminf.com/content/6/1/7>`_

..  _bidirectional:

Bidirectional Algorithm
~~~~~~~~~~~~~~~~~~~~~~~

The second algorithm included in the library is the bidirectional algorithm.
The goal is the same as in the classic approach -- to find
a path from the `source molecule` to the `target molecule`. However,
in this approach we use two trees (**A** and **B**) at once to get from M\ :sub:`S` to M\ :sub:`T`
(see :numref:`bidirectional-fig`). One is set up to search for a path from
M\ :sub:`S` to M\ :sub:`T` while the other searches in the opposite direction,
that is from M\ :sub:`T` to M\ :sub:`S`.

After every complete iteration of the original algorithm in both trees
(after the new structures are connected)
the `target molecule` in each tree is replaced by the closest
molecule to the current target from the opposite tree. This process is depicted
in :numref:`bidirectional-fig` where m\ :sub:`TA` is found to be
the closest molecule in tree **B** and, thus, becomes the new target of **A**.
Similarly, molecule m\ :sub:`TB` is the closest molecule to the current target
of **A** and, therefore, becomes the new target of **B**. This causes that
with each iteration the target of each tree becomes closer and closer.
The algorithm ends when any of the two trees finds its target
(the `connecting_molecule`).
This structure is guaranteed to exist in both trees and is used
to backtrack through them and generate the final path from M\ :sub:`S` to M\ :sub:`T`.

..      figure:: bidirectional.png
        :scale: 100%
        :name: bidirectional-fig

        Schematic depiction of the bidirectional algorithm. Two trees (**A** and **B**) are built
        at the same time with the `target <target molecule>`
        for each tree (m\ :sub:`TA` for **A** and m\ :sub:`TB` for **B**) being adjusted dynamically during runtime.

This algorithm should have some advantages over the classic one.
The biggest is probably the fact that the space
between M\ :sub:`S` and M\ :sub:`T` is explored in a more
uniform way. The classic algorithm can often converge
quickly to an area very close to M\ :sub:`T`, but
it usually takes a few more iterations before it finds the actual
structure of M\ :sub:`T`. This results in a disproportionate
number of molecules on the path being highly similar to M\ :sub:`T`.
This is not the case for the bidirectional search where
by the time the two trees get very close there
is already a lot of similar molecules between them
and it is much easier to find a connecting structure.

The bidirectional
approach might also be able to converge quicker,
because with each iteration the target gets closer and closer for each
tree and this significantly reduces the search space.

The following script shows how this algorithm can be used to generate
paths:

..  literalinclude:: ../../../src/python/molpher/examples/bidirectional.py
        :language: python
        :caption: Example implementation of a bidirectional pathfinding algorithm.
        :name: bidirectional-example
        :linenos:

The code above is almost identical to the example above (see :numref:`classic-example`).
It follows the same principles and the `Settings` class is also the same. The only difference
is that the run method is imported from `molpher.algorithms.bidirectional`
rather than `molpher.algorithms.classic`.

Antidecoys Algorithm
~~~~~~~~~~~~~~~~~~~~

.. todo:: Explain the algorithm and show pictures.

Stuff...