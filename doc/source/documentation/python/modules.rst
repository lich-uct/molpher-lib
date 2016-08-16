Python
======

This section documents the Python portion of the library. There are these
important packages under the main :mod:`molpher` package:

   #. :mod:`molpher.core`
         This is the package meant to be used by external scripts and applications.
         It contains a feature rich API with data structures and constructs
         that are more Python-like and easy to use.
   #. :mod:`molpher.algorithms`
         This package contains implementations of complete chemical space exploration algorithms.
         They are each located within their own subpackage with a :samp:`run` module
         that is used to set up and initiate the exploration.
   #. :mod:`molpher.swig_wrappers`
         Contents of this package are automatically generated from the C++ interface.
         Unlike the facilities under :mod:`molpher.core`, these are not intended for general use,
         but are used internally as a base for the constructs in :mod:`molpher.core`.
         However, they closely follow the C++ API so if you want to use the library
         from C++, you can get some idea on the data structures available from there
         as well.

There is also the `molpher.examples` which contains runnable scripts
that were discussed in the :ref:`tutorial`.

.. toctree::
   :maxdepth: 3

   molpher
