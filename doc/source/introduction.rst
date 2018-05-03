Introduction
============

This section attempts to introduce the uninformed to the main ideas behind :term:`molecular morphing`
and outline possible applications of this software. It also clarifies some terminology used in the documentation.

What is Molpher?
----------------

In early stages of a drug discovery endeavor a large library of compounds
is often screened for molecules that show some biological activity on a given target.
This is usually done in a laboratory setting using a technique called
:term:`high-throughput screening` (HTS).
HTS can be a very effective method, but is still quite limited in the number
of compounds that can be screened when we consider the incredible vastness of :term:`chemical space`
not to mention the costs associated with using and maintaining such sophisticated equipment
and actual physical library of compounds. Therefore, in practice, so called
:term:`virtual screening` (VS), is often employed prior to the actual physical experiments.
During VS, a virtual library of compounds is searched using computational methods,
which helps significantly in reducing the cost and time needed to find
new interesting structures of active molecules.

:term:`Molpher` is a tool developed for the purpose of designing a focused
virtual library of compounds that could be subsequently used in a virtual screening
campaign. Since Molpher creates compounds *de novo*, the generated library will very likely contain new chemistry
(unlike traditional VS libraries that are usually limited to known chemical space).
In Molpher, the goal of generating new compounds is achieved by connecting known active molecules
by a so called :term:`chemical space path`. Just like any path, a :term:`chemical space path`
has a start and a destination. In :term:`chemical space` this means
a :term:`source molecule` and a :term:`target molecule`. The path
is then just a sequence of structures that connect the two and is generated
by applying small structural changes to the source molecule.
These changes are called :term:`chemical operators` in the software and their
application is optimized to drive the search towards the target molecule.
The calculation ends when the target structure is generated and the path
is, thus, complete.

The main premise of Molpher is that the generated compounds, in particular those in the middle of the path, combine
the structural features of both the source and target. If the source and target are bioactive compounds, it is very likely
that some combination of their features will result in a bioactive compound as well.

..  note:: If you want to know more about how the original Molpher algorithm works, read `morphing-algorithm`
        section of the :doc:`tutorial <usage/tutorial>`.

What is Molpher-lib?
--------------------

:term:`Molpher` itself is a finished piece of software that efficiently implements
the basic molecular morphing algorithm described above, but has important limitations.
For example, there are many
possible implementations of the search and while many aspects of the
exploration can be parametrized, it is very hard to implement some
custom rules to prioritize certain compounds over others or
to provide customized operators,
which could be an especially interesting prospect when we think about synthesis
of compounds as well, for example.
In Molpher, there is also no way to maintain certain pharmacophore features,
which could dramatically increase
yield of the whole procedure when some mutual functional group arrangement is required.
Thus, came the motivation to develop a more flexible
and extensible solution that would give everyone an opportunity
to easily implement their own ideas tailored to their needs or just
integrate :term:`molecular morphing` with already existing solutions more easily.

To make it as straightforward as possible to get up to speed with the library,
we offer some guidance in the :doc:`next section <usage/index>`
on how to :doc:`install <usage/installation>` and :doc:`use <usage/tutorial>` the library as well as
an introduction into some new implementations of exploration :doc:`algorithms <usage/tutorial/algorithms>` in the library.