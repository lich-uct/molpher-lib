Introduction
============

This section attempts to introduce the uninformed to the main ideas behind the molecular morphing
method and the software. It also clarifies some terminology used in the documentation.

What is Molpher?
----------------

In early stages of any drug discovery endeavor a large library of compounds
often needs to be screened
for molecules that have the potential to become a new lead in the project. This
is usually done in a laboratory setting using a technique called
`high-throughput screening (HTS) <https://en.wikipedia.org/wiki/High-throughput_screening>`_.
HTS can be a very effective method, but is still very limiting in the number
of compounds that can be screened when we consider the incredible vastness of :term:`chemical space`.

:term:`Molpher` is a tool developed for the purpose
of finding  novel active compounds for synthesis/purchase and subsequent testing
in drug discovery projects.
It is capable of generating a library of potentially active compounds by
connecting known active molecules
with a :term:`chemical space path`. Just like any path, a :term:`chemical space path`
has a start and a destination. However, in :term:`chemical space` this means
a :term:`source molecule` and a :term:`target molecule`. The path
is then just a sequence of structures that connect the two (see :term:`chemical space path`
for more details).

The main motivation to generate such a path is
that some new interesting compounds might be discovered along the way.
The ones in the middle of the path are usually structurally different from both
the source and the target, but at the same time combine
their structural features, thus, some of them may form a good basis for
a new drug.

..  note:: If you want to know more about how the algorithm works, read `morphing-algorithm`
        section of the :doc:`tutorial <usage/tutorial>`.

What is Molpher-lib?
--------------------

:term:`Molpher` itself is aimed at users that might use it
in their day to day practice to generate focused virtual libraries of compounds.
However, there are many
possible implementations of the search and many aspects of the
exploration could be parametrized.
Thus, came the motivation to develop a more flexible
and extensible solution that would give everyone an opportunity
to easily implement their own ideas or just easily integrate :term:`molecular morphing`
in their existing solutions.

To make it as easy as possible to get up to speed with the library,
we have prepared a :doc:`tutorial <usage/tutorial>` in the :doc:`next section <usage/index>`
which contains information on how to install and use the library as well as
a few example implementations of exploration algorithms using the library.