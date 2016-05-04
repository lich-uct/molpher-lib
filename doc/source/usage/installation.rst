Installing the Library
======================

Below are two currently available methods of installation.

..  note:: As of yet, the library is only available for the Linux platform and will only run
        on 64-bit systems. However, other platforms will probably be supported in the future as well.

Installation with Anaconda
--------------------------

This is probably the easiest and preferred way to install the library and the associated Python package.
All you need to do is just either get the full `Anaconda <https://www.continuum.io/downloads>`_ distribution
or its lightweight variant, `Miniconda <http://conda.pydata.org/miniconda.html>`_.

Anaconda is essentially a Python distribution, package manager and virtual environment in one and makes setting up
a development environment for your project very easy.
After installing Anaconda/Miniconda all you need to do is to run the following
in the terminal:

..  code-block:: bash

    conda install -c lich molpher-lib

This will automatically download the latest version of the library and the corresponding Python package
and install in the currently active environment (for more information on environments and the
:command:`conda` command see `Conda Test Drive <http://conda.pydata.org/docs/test-drive.html>`_
in the official documentation).

Building and Installing from source
-----------------------------------

If you are a developer and want to build the most recent version from source,
you will need build everything from scratch. This section
describes how to do that on the Linux platform
and hints at possible problems one might encounter. If you want to build the library on
a different platform, this section might provide useful information as well.

Prerequisites
~~~~~~~~~~~~~

Before you start building, there are a few prerequisites that need to be satisfied on your system:

    - *build-essential* -- You will need this package to be able to build the software. It contains a compiler and other tools important for the build.
    - *swig* -- Used to generate the Python wrapping code (at least version 3.0 or later).
    - *bison* and *flex* -- Required by the rdkit toolkit.
    - *setuptools* and *pip* -- If you want to use *pip* to install the Python package, you will need this to be installed in your Python distribution.
    - *python{version}-dev* -- You will need this package in order to build the Python bindings for the particular Python *version*.

When those requirements are satisfied you can clone the source code from the repository:

..  code-block:: bash

    git clone https://{username}@bitbucket.org/Krwemil/molpher-lib.git # (replace {username} with your Bitbucket username)

and start building.

..  note:: The repository is currently private and you will need the authorization of the author to access it.

Building Molpher and Its Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, download the `Molpher Dependency Bundle
<https://drive.google.com/file/d/0B2rizkCQQcoybFdhOFExaVk5c0U/view?usp=sharing>`_
and extract the contents to a directory called :file:`deps/` of the repository root.
This bundle contains all the source code and binaries of the dependencies
that you will need to build and use the library on the Linux platform.

When it's done, just run:

..  code-block:: bash

    make CONF=Release_Linux_amd64

to build the release version of the library (no debugging symbols SWIG code is not generated) or:

..  code-block:: bash

    make CONF=Debug_SWIG

to build the development version with debug symbols and the SWIG code regenerated.

If you want to build the conda packages as well you can use the following
python script in the repository root:

..  code-block:: bash

    python build_conda.py

..  attention:: You will need the *jinja2* Python library to do that.

You can then install the packages by running:

..  code-block:: bash

    conda install --use-local molpher-lib

You can also install the Python package directly with *pip*:

..  code-block:: bash

    pip install .

in the repository root. This will install the Python package to the currently active Python distribution
and links the compiled extension against the compiled library in the :file:`lib/` directory generated
previously by :command:`make`.