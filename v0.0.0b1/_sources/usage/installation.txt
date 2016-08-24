Installation
============

Below are two currently available methods of installation.

..  note:: As of yet, the library is only available for the Linux platform and will only run
        on 64-bit systems. However, other platforms will probably be supported in the future as well.

..  warning:: The library has only been tested on a few Debian-based systems so far. Therefore, it might not work properly
    on other Linux distributions.

Installation with Anaconda
--------------------------

This is probably the easiest and preferred way to install the library and the associated Python package.
All you need to do is either get the full `Anaconda <https://www.continuum.io/downloads>`_ distribution
or its lightweight variant, `Miniconda <http://conda.pydata.org/miniconda.html>`_.

Anaconda/Miniconda is essentially a Python distribution, package manager and virtual environment in one and makes setting up
a development environment for your project very easy.
After installing Anaconda/Miniconda you need to run the following
to install the package:

..  code-block:: bash

    conda install -c lich molpher-lib

This will automatically install the latest non-development version of the library
in the currently active environment (for more information on environments and the
:command:`conda` command see `Conda Test Drive <http://conda.pydata.org/docs/test-drive.html>`_).

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

    - *build-essential* -- You will need this package in order to be able to build software on most Linux platforms. It contains a compiler and other tools important for the build.
    - *swig* -- It is used to generate the Python wrapping code (we are using version 3.0.10 at the moment).

        SWIG is only required if you made changes to the interface (header files under :file:`include/`) and want to configure cmake
        with the :command:`-DRUN_SWIG=ON` option (see the description of the *molpher_build_SWIG_Python* target in the section below).
        If this option is turned on, SWIG will be invoked by make upon build with the :command:`swig3.0`
        command so make sure the SWIG executable is available in the working environment.

    - *setuptools* -- This Python package is needed to install the associated Python code.
    - *python{version}-dev* -- You will need this package to build Python bindings for your Python *version*.

        If you get 'Missing Python.h' compiler errors, you probably do not have this package installed.

    - *dependencies* -- Molpher-lib depends on three third-party libraries:

        - *boost* (1.50.0)
        - *rdkit* (2014.03.1)
        - *tbb* (4.2)

        You can build the individual dependencies yourself and place them in the :file:`deps/` folder
        in the repository root. For each dependency, there should be a folder of the same name under :file:`deps/`
        (for example, the path to the *tbb* files would be :file:`deps/tbb/`). The cmake build script will automatically
        identify and prioritize dependencies in this directory. However, there is a build script (:file:`deps/build_deps.sh`)
        which can download and build the libraries automatically. Therefore, it should be sufficient to do:

        ..  code-block:: bash

            ./build_deps.sh --all

        You can also install the libraries on
        your system using a package manager or other means. In that case, cmake will automatically try to find them and
        link them during the build. Please, note that the build system is configured so that only static
        libraries of *boost* and *rdkit* are recognized. If you wish to link dynamically, you will need to modify
        the :file:`CMakeLists.txt` file accordingly.

Building the Library
~~~~~~~~~~~~~~~~~~~~

When the above requirements are met, you can start building. First, you need to initialize the cmake project
from a build directory:

..  code-block:: bash

    mkdir ${REPOSITORY_ROOT}/cmake_dir/ # create a subdirectory in the root of the repository
    cd ${REPOSITORY_ROOT}/cmake_dir/
    cmake .. # initialize the project

At this point, cmake checks the dependencies and generates the main makefile which defines three important targets:

    1. *molpher* -- This is the main target. It builds the C++ source code of the library and creates all necessary binary files. All targets mentioned below depend on this target.

    2. *molpher_install* -- This target will install the library in the given location.

        By default, this location is the :file:`dist/` folder in the repository root.
        This can be changed when the cmake project is initialized by setting
        `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/v3.5/variable/CMAKE_INSTALL_PREFIX.html>`_.
        Also by default, the TBB library is installed along with the molpher-lib files, because it
        is also a runtime dependency. If you do not want this behaviour you can instruct
        cmake to omit the installation of TBB with the :command:`-DINSTALL_TBB=OFF` option.

    3. *molpher_install_python* -- This will install the Python package and the C++ extensions using the generated SWIG wrapper code.

        By default, the primary Python distribution on the system is used. You can specify a different executable
        by specifying the path with :command:`-DPYTHON_EXECUTABLE`
        when you run cmake.

        If you want to update the SWIG wrapping code before this target is run, you can instruct cmake to do so with
        the :command:`-DRUN_SWIG=ON` option.

        When this target finishes, all required files should be in place and you should be able to
        import the *molpher* Python package located under :file:`src/python/`. Just add this directory to your :envvar:`PYTHONPATH` variable
        and it should work.

        ..  note:: You can also use the package created under :envvar:`CMAKE_INSTALL_PREFIX`. In this case, you need to have the following variables set up:

               ..  code-block:: bash

                    PYTHONPATH=$CMAKE_INSTALL_PREFIX/lib/pythonX.Y/site-packages # replace X.Y with your Python version
                    LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib

The build of a given configuration is initialized from the cmake project directory as:

..  code-block:: bash

    make $CONFIG # CONFIG is a configurations' name

Building Conda Packages
~~~~~~~~~~~~~~~~~~~~~~~

If you want to build your own conda packages, you can use a
python script located in :file:`conda/` subdirectory of the repository root:

..  code-block:: bash

    cd ${REPOSITORY_ROOT}/conda
    python build_conda.py

..  attention:: You will need `conda-build <https://github.com/conda/conda-build>`_ and the  *jinja2* Python library to do that.

You can install the built packages as follows:

..  code-block:: bash

    conda install --use-local molpher-lib

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

This documentation is generated using the :file:`build_docs.sh`
script under the :file:`doc/` directory. However, in order to successfully build the documentation
you will need a few packages in your Python environment:

..  literalinclude:: ../../../environment.yml
    :language: none
    :caption: The conda environment file used to build the documentation and test the library.

You can easily install all these packages like so:

..  code-block:: bash

    conda env create -n molpher-lib-docs -f environment.yml

Also, make sure you initialized the Python package (built the *molpher_install_python* target as described above)
prior to generating the documentation so that Sphinx can successfully import the *molpher* package.

To update the GitHub pages with the current version of the documentation, you can simply
run the :file:`build_docs.sh` with the ```--upload``` option:

..  code-block:: bash

    build_docs.sh --upload

..  note:: You will need write access to the repository and an SSH key attached to your account to be able to upload
        the documentation.
