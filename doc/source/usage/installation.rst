Installation
============

In this section, we discuss various methods of installation. There are also instructions on how to build the library
from source and other comments useful to library developers.

..  note:: As of yet, the library is only available for the Linux platform and will only run
        on 64-bit systems. However, support for other platforms is on the roadmap as well.

..  warning:: So far, the library is only tested on a limited number of Debian-based systems. Therefore, it might not work properly
    on some other Linux distributions.

Installation with Anaconda
--------------------------

If you want to use the library from Python environment, this is probably the best way for you to install.
All you need to do is either get the full `Anaconda <https://www.continuum.io/downloads>`_ distribution
or its lightweight variant, `Miniconda <http://conda.pydata.org/miniconda.html>`_.

Anaconda/Miniconda is essentially a Python distribution, package manager and virtual environment in one and makes setting up
a development environment for your projects very easy.
After installing Anaconda/Miniconda you need to run the following
to install the package in your environment:

..  code-block:: bash

    conda install -c rdkit -c lich molpher-lib

This will automatically install the latest non-development version of the library
in your currently active environment.
If you are interested in the development snapshots, you can specify the :code:`dev` label while you install:

..  code-block:: bash

    conda install -c rdkit -c lich/label/dev molpher-lib

..  attention:: The binaries in the :code:`molpher-lib` package are compiled against
    the newest version of RDKit (:code:`2018.03.1` at the time of writing),
    which has a different ABI than the previous versions.
    Therefore, this and the following versions of Molpher-lib will not work with older versions of RDKit anymore.

For more information on environments and the
:command:`conda` command see `Conda Test Drive <http://conda.pydata.org/docs/test-drive.html>`_.

You can test if the library works as it should by running the Python test suite:

..  code-block:: python

    from molpher.tests import run

    run()

Building and Installing from Source -- Linux
--------------------------------------------

If you want to build the source code yourself,
you will have to use :command:`cmake`. This section
describes how to do that on the Linux platform. If you want to build the library on
a different platform, this section might provide useful information as well. However,
you might need to make some adjustments to the underlying :file:`CMakeLists.txt` file
and build and link the dependencies yourself.

..  note:: If you manage to successfully build the library on a platform other than Linux,
        please, consider making a pull request and contributing your code to the `official GitHub
        repository <https://github.com/lich-uct/molpher-lib.git>`_.

You can get the source code with the following command:

..  code-block:: bash

    git clone https://github.com/lich-uct/molpher-lib.git

Prerequisites
~~~~~~~~~~~~~

Before you start building, there are a few requirements that need to be satisfied:

    - *cmake* -- This tool will generate the Makefiles for the project and is used to configure the build. Cmake version 3.9 and higher is supported.
    - *build-essential* -- You will need this package in order to be able to build software on most Linux platforms. It contains a compiler and other tools important for the build.

        Note that you might encounter some problems if you decide to use an ancient compiler with poor
        support for the newer C++ standards (C++11 onwards). For example, g++ 5.4 and newer should be OK,
        but even slightly older compilers could work with no problems.

    - *swig* -- It is used to generate the Python wrapping code (we are using version 3.0.12 at the moment).

        SWIG is only required if you made changes to the binary interface (header files under :file:`include/`) and want to configure cmake
        with the :command:`-DRUN_SWIG=ON` option (see the description of the *molpher_build_SWIG_Python* target in the section below).
        If this option is turned on, SWIG will be invoked by make upon build with the :command:`swig3.0`
        command so make sure the SWIG executable is available in the working environment.

    - *setuptools* -- This Python package is needed to build and install the Molpher-lib Python package.
    - *python{version}-dev* -- You will need this package to build Python bindings for your Python *version*.

        If you get 'Missing Python.h' compiler errors, you probably do not have this package installed.

    - *dependencies* -- Molpher-lib depends on three third-party libraries:

        - *tbb* (most versions should work fine, we generally build against 2018 Update 3)
        - *boost* (most versions should work fine, we generally build against 1.65)
        - *rdkit* (2018.03.1 and newer)
        - *numpy* (RDKit dependency in Python, not required if you will be using the C++ interface only)

        There is a bash script (:file:`deps/build_deps.sh`)
        which can download and build the dependencies automatically with
        the required options. It should be sufficient to just run:

        ..  code-block:: bash

            ./build_deps.sh --all

        If you want the dependencies yourself, you should install them in the :file:`deps/` folder
        in the repository root. For each dependency, there should be a folder of the same name under :file:`deps/`
        (for example, the path to the *tbb* files would be :file:`deps/tbb/`). The CMakeLists.txt is configured to automatically
        identify and prioritize dependencies in this directory.

        You can also leverage the libraries already installed on your system. In that case, cmake should automatically find them and
        link them during the build. The :file:`CMakeLists.txt` file is configured to link against dynamic versions of
        all libraries so make sure you have those installed.

Building the Library
~~~~~~~~~~~~~~~~~~~~

When the above requirements are met, you can start building. First, you need to create a build directory
and initialize the cmake project:

..  code-block:: bash

    mkdir ${REPOSITORY_ROOT}/cmake-build/ # create a subdirectory in the root of the repository
    cd ${REPOSITORY_ROOT}/cmake-build/
    cmake ..

This is the simplest configuration with default options, but most of the time we will probably
require more customization. The cmake configuration file recognizes a few options.
For example, the following will force debug mode and Python 3 during build:

..  code-block:: bash

    cmake .. -DCMAKE_BUILD_TYPE=Debug -DPYTHON_EXECUTABLE=python3

If you want to recreate the Python wrapping code during build, you should
add :command:`-DRUN_SWIG=ON`. Remember, that you need to have SWIG installed in a standard
location for this to work. Alternatively, you can add swig to your :envvar:`PATH` or use
:command:`-DSWIG_EXECUTABLE=/path/to/swig` to tell cmake where to look for it.

When the makefile is created, you can use :command:`make` to build various targets:

..  code-block:: bash

    make $CONFIG # CONFIG is a configurations' name

There are three important targets:

    1. *molpher* -- Builds the binaries for the C++ part of Molpher-lib.

    2. *molpher_install* -- This target will install the library in a given location.

        By default, this location is the :file:`dist/` folder in the repository root.
        This can be changed when the cmake project is initialized by setting
        `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/v3.9/variable/CMAKE_INSTALL_PREFIX.html>`_.
        By default, the required dependency libraries are not installed. If you just want to install them,
        you can configure cmake to do so by setting the following: :command:`-DINSTALL_TBB=ON -DINSTALL_Boost=ON -DINSTALL_RDKit=ON`.

    3. *molpher_install_python* -- This builds the C++ Python extension and installs the Python package into :envvar:`CMAKE_INSTALL_PREFIX`.

        By default, the primary Python distribution on the system is used. You can specify a different executable
        by with :command:`-DPYTHON_EXECUTABLE`.

        If you want to update the SWIG wrapping code before this target is run, you can instruct cmake to do so with
        the :command:`-DRUN_SWIG=ON` option. Do not forget to specify the swig path with :command:`-DSWIG_EXECUTABLE` if it is installed in a non-standard location.

        When this target finishes, all required files should be in place and you should be able to
        import the *molpher* Python package, provided that your :envvar:`PYTHONPATH` and :envvar:`LD_LIBRARY_PATH` are set
        appropriately. Here is an example of how these variables can be set if standard locations are used:

        ..  code-block:: bash

            export CMAKE_INSTALL_PREFIX="${REPOSITORY_ROOT}/dist"
            export DEPS_DIR=${CMAKE_INSTALL_PREFIX}/../deps
            export PYTHONPATH=${DEPS_DIR}/rdkit/:${CMAKE_INSTALL_PREFIX}/lib/python3.5/site-packages
            export LD_LIBRARY_PATH=${DEPS_DIR}/tbb/lib/intel64/gcc4.7:${DEPS_DIR}/rdkit/lib/:${DEPS_DIR}/boost/stage/lib:${CMAKE_INSTALL_PREFIX}/lib

        You should then be able to successfully run the Python test suite:

        ..  code-block:: python

            from molpher.tests import run

            run()

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

This documentation is generated using the :file:`build_docs.sh`
script under the :file:`doc` directory. However, you will need a few Python packages
in order to successfully build it. Molpher-lib source code contains a conda :download:`environment file <../../../environment.yml>`
which defines these requirements. You can install this environment like so:

..  code-block:: bash

    conda env create -f environment.yml

The resulting environment will be called *molpher-lib* and you can activate it while setting
important Molpher-lib variables by sourcing the :file:`source_2_activate`
file in :envvar:`REPOSITORY_ROOT`:

..  code-block:: bash

    . source_2_activate

This will not only allow you to build
the documentation, but also run the code in the associated Jupyter notebooks.

Once your environment is activated, you can build the Python wrappers and
generate the documentation:

..  code-block:: bash

    python setup.py build_ext --inplace
    cd doc
    ./build_docs.sh

To update the GitHub pages, it is possible to
run the :file:`build_docs.sh` script with the ```--upload``` option:

..  code-block:: bash

    build_docs.sh --upload

..  note:: You will need write access to the repository and an SSH key attached to your account to be able to do this.

Building Conda Packages
~~~~~~~~~~~~~~~~~~~~~~~

If you want to build your own conda packages, you can use a
python script located in :file:`conda` subdirectory of the repository root:

..  code-block:: bash

    cd ${REPOSITORY_ROOT}/conda
    python build.py

..  attention:: You will need `conda-build <https://github.com/conda/conda-build>`_ and the  *jinja2* Python library to do that.
    These are both part of the *molpher-lib* conda environment we introduced before. It is enough to just do :code:`conda activate molpher-lib`.

The built packages will be located at :file:`/tmp/conda-bld`.
