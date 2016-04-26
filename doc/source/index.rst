
Molpher-lib: Python and C++ API for Molpher
===========================================

What you have just stumbled upon is the documentation of *Molpher-lib*,
a computer library for exploration of chemical space with :term:`molecular morphing`.
It is essentially a fork of :term:`Molpher`
-- a free open source software developed for the same purpose.

The original :term:`Molpher` is mostly aimed at users. Therefore, it contains a rich
graphical user interface with visualizations of chemical space and many facilities
to directly affect the exploration process.
However, it is not very friendly
when it comes to developers
who usually prefer access to the facilities of a computer program via an API
or would even like to alter the implementation so that it suits their needs.
We have become more and more aware of these shortcomings during our research
where we would often like to implement new ways of filtering the generated :term:`morphs <morph>`
or evaluating their 'fitness', but the monolithic implementation in C++ would make
it difficult and time-consuming.

Therefore, the goal of this project is to make :term:`Molpher` more accessible to
programmers, developers and scientists alike
by providing a flexible and extensible programming interface including
a Python package that exposes most of the features of *Molpher-lib* to the
Python programming language as well.

..  warning:: Please, note that the project is in early stage of development
        and that a lot of the features mentioned in this documentation might
        change in the future and that there might still be a considerable amount
        of bugs and other issues.

        .. todo:: Add link to the issue tracker

Table of Contents
=================

..  toctree::
    :name: maintoc
    :numbered:
    :maxdepth: 5

    self
    Introduction <introduction>
    usage/index
    documentation/index

Indices and Search
==================
* :ref:`modindex`
* :ref:`genindex`
* :ref:`glossary`

TODO list
=========

This documentation should contain most of the information an average user might need,
but there are still places where it needs to be clarified or completed. Here is
a list of items that need further revision:

.. todolist::
