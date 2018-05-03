..  _morphing-basics:

Molecular Morphing Basics
=========================

This section showcases the most basic features and building blocks
needed to successfully start modifying compounds in an automated fashion
with the library. It also hints at possible applications
that can go beyond molecular morphing.

Not all of the features and settings will be presented and we encourage
you to take a look at the `source-code-docs` if you want to know more about the implementation itself
or more advanced topics.

The Python examples in this tutorial are taken from a sample Jupyter Notebook which
can be downloaded :download:`here <../../../notebooks/basics.ipynb>` and the rendered HTML version
can be viewed `here </_static/basics.html>`_.

If you run into an issue with the examples or have an idea or a suggestion,
consider submitting to the `issue tracker <https://github.com/lich-uct/molpher-lib/issues>`_.

.. _morphing-algorithm:

The Molecular Morphing Concept
------------------------------

Let us first shortly describe how :term:`molecular morphing` works. This should help anyone who has never heard
about the method to understand it and get an idea on what this library can do for them.

In order to find a path that connects the :term:`source molecule` and the :term:`target molecule`
(see `../../introduction` if you are not familiar with these terms),
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