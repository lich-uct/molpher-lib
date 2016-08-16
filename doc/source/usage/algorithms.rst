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

This is the algorithm published along with the :term:`Molpher`
program [1]_.

.. todo:: Explain the algorithm and show pictures.

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

The second example we have here is a little bit more elaborate,
but implements a very simple idea. Instead of one exploration tree,
we build two trees that each search for a path to the closest molecule in the other:

.. todo:: Explain the algorithm and show pictures.

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

Antidecoys Algorithm
~~~~~~~~~~~~~~~~~~~~

.. todo:: Explain the algorithm and show pictures.

Stuff...