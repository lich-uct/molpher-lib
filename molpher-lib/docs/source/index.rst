.. molpher-lib documentation master file, created by
   sphinx-quickstart on Thu Nov 26 10:28:29 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Molpher-lib: Python and C++ API for Molpher
===========================================

What you have just stumbled upon is the documentation of *Molpher-lib*.
*Molpher-lib* is essentially a C++ API derived from :term:`Molpher`
-- a free open source software aimed at exploration of `chemical space`.
It shouldn't be regarded as a full API to :term:`Molpher`, but as a
standalone project that just builds upon its features.
However, it is not out of the question that both of these projects
will get more tightly integrated in the future.

The goal of this project is to make :term:`Molpher` more accessible to developers
by providing a flexible and extensible programming interface.
It also tries to make the functionality available for wider audience by
including a Python package that exposes the API to Python
using the `Simplified Wrapper and Interface Generator (SWIG) <http://swig.org/>`_.

..  warning:: Please, note that the project is in a very early stage of developement
        and that a lot of features mentioned in this documentation might, and probably
        will, change in the future. Additionally, be wary of any use in production, because
        bugs and performance issues will most probably occur.

Table of Contents
=================

..  toctree::
    :name: maintoc
    :numbered:
    :maxdepth: 2

    self
    Introduction <introduction>
    usage/index
    documentation/index

Indices and Search
==================
* :ref:`search`
* :ref:`modindex`
* :ref:`genindex`
* :ref:`glossary`

TODO list
=========

This documentation should contain most of the information an average user might need,
but there are still places where it needs to be clarified or completed. Here is
a list of items that need further revision:

.. todolist::
