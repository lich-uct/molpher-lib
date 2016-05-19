Installing the Library
======================

Below are two currently available methods of installation.

..  note:: As of yet, the library is only available for the Linux platform and will only run
        on 64-bit systems. However, other platforms will probably be supported in the future as well.

..  warning:: The library has only been tested on a few Debian-based systems so far. Therefore, it might not work on
    different Linux distributions.

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
you will need to build everything from scratch. This section
describes how to do that on the Linux platform
and hints at possible problems one might encounter. If you want to build the library on
a different platform, this section might provide useful information as well.

You can get the source code with the following git command:

..  code-block:: bash

    git clone https://github.com/lich-uct/molpher-lib.git

Prerequisites
~~~~~~~~~~~~~

Before you start building, there are a few requirements that need to be satisfied:

    - *build-essential* -- You will need this package to be able to build the software. It contains a compiler and other tools important for the build.
    - *swig* -- It is used to generate the Python wrapping code (we are using version 3.0.8 at the moment).

        Only required if you made changes to the interface (header files under :file:`include/`) and need to compile the code
        using the :command:`Debug_SWIG` configuration. SWIG is invoked by make upon build with the :command:`swig3.0`
        command so make sure it resides somewhere on your :makevar:`PATH`.

    - *setuptools* and *pip* -- If you want to use *pip* or *easy_install* to install the Python package, you will need this to be included in your Python distribution.
    - *python{version}-dev* -- You will need this package in order to build the Python bindings for the particular Python *version*.

        If you get 'Missing Python.h' errors, you probably do not have this package installed.

    - *dependencies* -- Molpher-lib depends on three third-party libraries:

        - *boost* (1.49.0)
        - *rdkit* (2014.03.1)
        - *tbb* (4.2)

        You can build the individual dependencies yourself and place them in the :file:`deps/` folder
        in the repository root. For each dependency, there should be a folder of the same name under :file:`deps/`
        (for example, the path to the *tbb* files would be :file:`deps/tbb/`).

        Alternatively, you can download the `Molpher Dependency Bundle
        <https://drive.google.com/file/d/0B2rizkCQQcoybFdhOFExaVk5c0U/view?usp=sharing>`_
        and extract the contents to the :file:`deps/` directory.
        This bundle contains pre-built dependencies for the Linux platform
        and should work out of the box on Debian-based systems.

Building the Library
~~~~~~~~~~~~~~~~~~~~

When the above requirements are met, you can start building. There are three configuration options:

    1. Debug -- Builds the library with debugging symbols, but without updating the SWIG wrapping code.

        If you did not make any changes to the header files under :file:`include/`,
        and do not require changes to the wrapping code, you can use this option.

    2. Debug_SWIG -- This builds the library and updates the SWIG wrapping code.

        Use it if you made
        changes to the header files under :file:`include/` or want to update the SWIG wrappers.

    3. Release_Linux_amd64 -- This builds the library for general use.

        Binaries produced with this
        configuration are intended for distribution.

You run the build process with:

..  code-block:: bash

    make CONF=$CONFIG_NAME # CONFIG_NAME is the name of one of the configurations above

This will compile the shared object file and put it in the :file:`dist/lib/`
directory along with the required *tbb* libraries that are dynamically
linked during runtime.

Building the Python Wrappers and the Conda Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The above procedure only builds the library itself. If you want to compile the Python bindings,
you can do so by invoking:

    ..  code-block:: bash

        python setup.py build_ext

in the repository root.

You can also install the Python package with, for example:

    ..  code-block:: bash

        pip install .

This will compile the Python extension and install the Python package to the currently active Python environment.

    .. note:: The extension is linked against the compiled library in the :file:`lib/` directory generated
        previously by :command:`make`.

If you want to build the conda packages as well, you can use the following
python script in the repository root:

..  code-block:: bash

    python build_conda.py

..  attention:: You will need `conda-build <https://github.com/conda/conda-build>`_ and the  *jinja2* Python library to do that.

You can install the built packages by running:

..  code-block:: bash

    conda install --use-local molpher-lib

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

The library has an extensive documentation that is generated using the :file:`build_docs.sh`
script under the :file:`doc/` directory. In order to successfully build the documentation,
you will need a few packages in your Python environment:

..  literalinclude:: ../../../doc/requirements.txt
    :language: none
    :caption: The requirements file for the environment used to build the documentation.

You can easily install the packages with *pip*:

..  code-block:: bash

    pip install -r requirements.txt

To update the GitHub pages with the current version of the documentation, you can simply
run the :file:`build_docs.sh` with the ```--upload``` option:

..  code-block:: bash

    build_docs.sh --upload

..  note:: You will need write access to the repository and an SSH key attached to your account to do that.
