Installation
============

In this section, we discuss various methods of installation. You will also find instructions on how to build the library
from source and create conda packages.

..  note:: As of yet, the library is only available for Linux. However, support for other platforms is
        on the roadmap as well.
        If you are interested in a build for a different platform, feel free to express your
        wish on the `issue tracker <https://github.com/lich-uct/molpher-lib/issues>`_.

Installation with Anaconda
--------------------------

This is probably the best way for you to install if you want easy access to both the C++ and Python parts of the library.
All you need to do is either get the full `Anaconda <https://www.continuum.io/downloads>`_ distribution
or its lightweight variant, `Miniconda <http://conda.pydata.org/miniconda.html>`_.

Anaconda/Miniconda is essentially a Python distribution, package manager and virtual environment in one and makes setting up
a development environment for your projects more straightforward.
After installing Anaconda/Miniconda you need to run the following command
to install the :code:`molpher-lib` package in your environment:

..  code-block:: bash

    conda install -c conda-forge -c lich molpher-lib

This should be all you need. The command will automatically install the
latest stable version of the library in the currently active conda environment.
If you want the newest features, you can also try the development snapshots
by specifying the :code:`dev` label while you install:

..  code-block:: bash

    conda install -c conda-forge -c lich/label/dev molpher-lib

Development snapshots have the newest functionality, but may be more buggy.

You can test if the library works as it should by running the Molpher-lib test suite in your Python shell:

..  code-block:: python

    from molpher.tests import run

    run()

Building and Installing from Source
-----------------------------------

If you are a developer and you want to modify and build the source code yourself,
you will have to use :command:`cmake` to build the library and the Python bindings. This section
describes how to do that on a Linux machine. If you want to build the library on
a different platform, this section will provide useful information as well.
You will need to make some adjustments to the underlying :file:`CMakeLists.txt` to support your platform, however.

..  note:: If you manage to successfully build the library on a platform other than Linux,
        please, consider making a pull request and contributing your code to the `official GitHub
        repository <https://github.com/lich-uct/molpher-lib.git>`_. Also consider giving `this issue <https://github.com/lich-uct/molpher-lib/issues/7>`_
        a bump up if you are interested in Windows build support.

Prerequisites
~~~~~~~~~~~~~

In order to build the library, some tools and libraries need to present.
You can learn about their role in detail here, but the actual setup can be done without much effort inside a conda environment so feel
free to skip to the :ref:`next section <conda-build-env>`

Overview of the tools involved (all can be obtained from Anaconda):

    - *git* -- Used to check out the source code from the repository.
    - *cmake* -- Generate the Makefiles for the project and configure the build (searches for dependencies, makes sure libraries are correctly linked, etc.). You should use version 3.9 and higher.
    - *build-essential* -- On most Linux distros, this package provides a compiler and other tools important for the build.
    - *swig* -- It is used to generate the Python wrapping code (we are using version 3.0 at the moment).

        SWIG is only required if you made changes to the binary interface (header files under :file:`include/`) and want to configure cmake
        with the :command:`-DRUN_SWIG=ON` option to build the Python bindings as well
        (see the description of the *molpher_install_python* target in the section below).
        If this option is turned on, SWIG will be invoked by make upon build with the :command:`swig`
        command so make sure the SWIG executable is available in the working environment or that you have
        correctly specified the :command:`-DCONDA_PREFIX` directive or :command:`CONDA_PREFIX` environment variable.

    - *setuptools* -- Needed to build and install the Molpher-lib Python package.
    - *python{version}-dev* -- You will need this package to build Python bindings for various Python versions.

        If you get 'Missing Python.h' compiler errors, you probably do not have this package installed.

    - *dependencies* -- Molpher-lib depends on three third-party libraries:

        - *tbb* (most versions should work fine up to 2020, versions starting 2021 lack some of the older interfaces and are currently not supported).
        - *boost* (newest versions should work fine, i.e. 1.74.0)
        - *rdkit* (there were some changes to the names of RDKit libraries and only versions newer than :command:`2019.03.4` are supported
        - *numpy* (RDKit dependency in Python, not required if you are only building the C++ code)

        You can leverage the libraries already installed on your system. In that case, :command:`cmake` should automatically find them
        on your path and link them during the build. The :file:`CMakeLists.txt` file is configured to link against dynamic versions of
        all libraries so make sure you have those installed. You can also obtain these dependencies in the conda build environment (see :ref:`*Using Conda to
        Manage the Build Environment* <conda-build-env>`).

.. _conda-build-env:

Using Conda to Manage the Build Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most straightforward way to setup a build environment for Molpher-lib is through Anaconda. You can start with the
:code:`environment.yml` file in the root directory of the project repository
(available `here <https://github.com/lich-uct/molpher-lib/blob/master/environment.yml>`_) and create a new
environment from it:

..  code-block:: bash

    conda env create -f environment.yml # the environment will be called molpher-lib-build by default

The resulting build environment will have all libraries and tools needed to build and install Molpher-lib.
You can then point :command:`cmake` to
the environment prefix with :command:`-DCONDA_PREFIX=/path/to/your/env/` or
just activate the environment, which will set the :command:`CONDA_PREFIX` variable automatically
and also provide you with a compatible version of :command:`cmake`.

Building the Library
~~~~~~~~~~~~~~~~~~~~

When your environment is set, you can start building. If you have not done so already, you need to
check out the code and create a build directory

..  code-block:: bash

    git clone https://github.com/lich-uct/molpher-lib.git
    REPOSITORY_ROOT="`pwd`/molpher-lib"
    mkdir ${REPOSITORY_ROOT}/cmake-build/

Then you can initialize the cmake project:

..  code-block:: bash

    # using conda for the build environment here
    conda activate molpher-lib-build

    cd ${REPOSITORY_ROOT}/cmake-build/
    cmake ..

This is the simplest configuration with default options, but sometimes you may
require more settings. The Molpher-lib configuration recognizes a few options
(take a look at :file:`CMakeLists.txt` and you will be able to see them).
For example, the following will force debug mode and Python 3 during build,
which is useful when you want to step through the code as you run it:

..  code-block:: bash

    cmake .. -DCMAKE_BUILD_TYPE=Debug

If you want to recreate the Python wrapping code during build (you changed the binary interface under
${REPOSITORY_ROOT}/include), you should
also add :command:`-DRUN_SWIG=ON`.

When the makefiles are created, you can use :command:`make` for different build and installation scenarios:

..  code-block:: bash

    make $CONFIG # $CONFIG is a configuration scenario

In the code above you can specify the following as :command:`$CONFIG`:

    1. *molpher* -- Builds only the C++ code.

    2. *molpher_install* -- Will install the binaries of the library in the given location.

        By default, this location is the :file:`dist/` folder in the repository root.
        This can be changed when the cmake project is initialized by setting
        `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/v3.9/variable/CMAKE_INSTALL_PREFIX.html>`_.
        By default, the required dependency libraries are not installed. If you want to install them with the library,
        you can configure cmake to do so by setting the following: :command:`-DINSTALL_TBB=ON -DINSTALL_Boost=ON -DINSTALL_RDKit=ON`.

    3. *molpher_install_python* -- This builds the C++ Python extension and installs the Python package into :envvar:`CMAKE_INSTALL_PREFIX`.

        By default, the Python distribution in the active environment is used. You can specify a different executable
        with :command:`-DPYTHON_EXECUTABLE`.

        If you want to update the SWIG wrapping code before this target is run, you can instruct cmake to do so with
        the :command:`-DRUN_SWIG=ON` option. Do not forget to specify the swig path with :command:`-DSWIG_EXECUTABLE` if it is installed in a non-standard location.

        When this target finishes, all required files should be in place and you should be able to
        import the *molpher* Python package, provided that your :envvar:`PYTHONPATH` and :envvar:`LD_LIBRARY_PATH` are set
        appropriately. Here is an example of how these variables can be set if a conda build environment is used to supply dependencies:

        ..  code-block:: bash

            # inside the repository root directory
            conda activate molpher-lib-build
            export PYTHONPATH="`pwd`/src/python/"
            export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/:`pwd`/dist/lib/"

            python # molpher package should now be importable

        You can test the package by running the Python test suite:

        ..  code-block:: python

            from molpher.tests import run

            run()

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

This documentation was generated with the :file:`build_docs.sh`
script under the :file:`doc` directory. The packages required
are also available in the :download:`environment file <../../../environment.yml>`
so after the built with *molpher_install_python* configuration,
you can also compile the documentation:

..  code-block:: bash

    cd doc
    ./build_docs.sh

To update the GitHub pages, it is possible to
run the :file:`build_docs.sh` script with the ```--upload``` option:

..  code-block:: bash

    build_docs.sh --upload

..  note:: You will need write access to the repository to be able to do this.

Building Conda Packages
~~~~~~~~~~~~~~~~~~~~~~~

If you want to package the build as a conda package, you can use a
python script located in the :file:`conda` subdirectory of the repository root:

..  code-block:: bash

    python build.py

..  note:: You will need `conda-build <https://github.com/conda/conda-build>`_ and the  *jinja2* Python library,
    but both are also included in the default environment file.

The packages are built for each version of Python separately and will be located at :file:`${REPOSITORY_ROOT}/conda-build/conda`.
From there, they can be uploaded to Anaconda Cloud.
