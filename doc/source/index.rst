
Molpher-lib: Systematic Chemical Space Exploration with Molecular Morphing
==========================================================================

This library is essentially a fork of the *de novo* drug design framework :term:`Molpher`. Molpher is
a program developed to tackle systematic exploration of chemical space
using an intuitive atom-based method called :term:`molecular morphing`
to sample new and potentially useful chemical structures. The task of Molpher is to bring the idea of molecular
morphing to non-technical users. Therefore, it contains a rich
graphical user interface with visualizations of chemical space and other GUI features that help
to analyse and oversee the exploration.

Over the years, however, as we learnt more about this concept
and its potential, the need to adjust certain parts of
the exploration algorithm for different purposes became more imminent and Molpher-lib was born. It wraps the core
idea of Molpher in a user-friendly programming interface (also available from Python)
that can be used to drive the
exploration in a more nuanced way or to implement entirely new use cases.

If you wish to jump right in, you can take a look at the :doc:`tutorial <usage/tutorial>`. It showcases some
of the typical tasks that the library was designed for and should get any new users acquainted with it quickly. If you wish
to know more about molecular morphing itself and how this library can be useful in *de novo* drug design,
you may continue reading the next :doc:`section <introduction>`. To understand the concept of molecular morphing more in depth,
you can check out the `original paper <https://dx.doi.org/10.1186/1758-2946-6-7>`_ on the Molpher software.

..  warning:: Please, note that the project is in early stage of development
        and so a lot of features mentioned in this documentation might
        change in the future and there might still be a considerable amount
        of bugs and other issues.

..  note:: If you find a bug or have a point to make
        regarding missing features or other things, do not hesitate to post
        on the project's `issue tracker <https://github.com/lich-uct/molpher-lib/issues>`_.

Table of Contents
=================

..  toctree::
    :name: maintoc
    :maxdepth: 5

    self
    Molecular Morphing Introduction <introduction>
    usage/index
    documentation/index
    glossary

Indices and Search
==================
* :ref:`modindex`
* :ref:`genindex`
* :ref:`glossary`

.. todolist::
