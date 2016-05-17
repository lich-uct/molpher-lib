
Molpher-lib: a C++/Python Library For Systematic Chemical Space Exploration
===========================================================================

This library is a fork of :term:`Molpher`, a program for exploration of
chemical space with :term:`molecular morphing`.
Molpher mostly targets non-technical users. Therefore, it contains a rich
graphical user interface with visualizations of chemical space and supports
various platforms.

However, Molpher is not very friendly
when it comes to developers
who usually require access to a computer program via an API
or would like to easily alter the implementation so that it suits better to
their needs. We have become more and more aware of the shortcomings of Molpher during our research
where we would often like to implement new ways of filtering the generated :term:`morphs <morph>`
or evaluating their 'fitness', but there was no easy way to do that,
because the underlying C++ implementation is not very flexible and user-friendly.

The goal of this project is to make molecular morphing more accessible to
programmers, developers and scientists alike
by providing a flexible and extensible programming interface that can also be
used from a high-level programming language such as Python.

..  warning:: Please, note that the project is in early stage of development
        and so a lot of features mentioned in this documentation might
        change in the future and there might still be a considerable amount
        of bugs and other issues.

        .. todo:: Add link to the issue tracker

Table of Contents
=================

..  toctree::
    :name: maintoc
    :maxdepth: 5

    self
    Introduction <introduction>
    usage/index
    documentation/index
    glossary

Indices and Search
==================
* :ref:`modindex`
* :ref:`genindex`
* :ref:`glossary`

.. todolist::
