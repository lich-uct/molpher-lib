Installing the Library
======================

Below are two currently available methods of installation.

..  note:: As of yet, the library is only available for the Linux platform and will only run
        on 64-bit systems. However, other platforms will probably be supported in the future as well.

..  warning:: The library has only been tested on a few Debian-based systems so far. Therefore, it might not work properly
    on different Linux distributions.

Installation with Anaconda
--------------------------

This is probably the easiest and preferred way to install the library and the associated Python package.
All you need to do is just either get the full `Anaconda <https://www.continuum.io/downloads>`_ distribution
or its lightweight variant, `Miniconda <http://conda.pydata.org/miniconda.html>`_.

Anaconda/Miniconda is essentially a Python distribution, package manager and virtual environment in one and makes setting up
a development environment for your project very easy.
After installing Anaconda/Miniconda all you need to do is to run the following
in the terminal:

..  code-block:: bash

    conda install -c lich molpher-lib

This will automatically download the latest version of the library and the corresponding Python package
and install in the currently active environment (for more information on environments and the
:command:`conda` command see `Conda Test Drive <http://conda.pydata.org/docs/test-drive.html>`_
in the official documentation).

Building and Installing from Source -- Linux
--------------------------------------------

If you want to build the source code yourself,
you can accomplish that quite easily using :command:`cmake`. This section
describes how to do that on the Linux platform. If you want to build the library on
a different platform, this section might provide useful information as well. However,
you might need to make some adjustments to the underlying :file:`CMakeLists.txt` file
and build and link the dependencies yourself.

..  note:: If you manage to successfully build the library on a platform other than Linux,
        please, consider making a pull request and contributing your code to the `official GitHub
        repository <https://github.com/lich-uct/molpher-lib.git>`_.

You can get the source code with the following git command:

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

    - *setuptools* and *pip* -- If you want to use *pip* or *easy_install* to install the Python package when the build and installation of the library are complete, you will need these to be installed with Python.
    - *python{version}-dev* -- You will need this package in order to build the Python bindings for the particular Python *version*.

        If you get 'Missing Python.h' compiler errors, you probably do not have this package installed.

    - *dependencies* -- Molpher-lib depends on three third-party libraries:

        - *boost* (1.49.0)
        - *rdkit* (2014.03.1)
        - *tbb* (4.2)

        You can build the individual dependencies yourself and place them in the :file:`deps/` folder
        in the repository root. For each dependency, there should be a folder of the same name under :file:`deps/`
        (for example, the path to the *tbb* files would be :file:`deps/tbb/`). The cmake build script will automatically
        identify and prioritize dependencies in this directory.

        You can also install the libraries on
        your system using a package manager or other means. In that case, cmake will automatically try to find them on
        your system and link them during the build. Please, note that the build system is configured so that only static
        libraries of *boost* and *rdkit* are recognized. If you wish to link dynamically, you will need to modify
        the :file:`CMakeLists.txt` file.

        If you are building on a 64-bit Linux machine, you can download the `Molpher Dependency Bundle
        <https://drive.google.com/file/d/0B2rizkCQQcoybFdhOFExaVk5c0U/view?usp=sharing>`_.
        This bundle contains pre-built dependencies for this platform
        and should work out of the box on Debian-based systems. You just need to extract
        the contents of the downloaded archive into the :file:`deps/` directory in the repository root.

Building the Library
~~~~~~~~~~~~~~~~~~~~

When the above requirements are met, you can start building. First, you need to initialize the cmake project
from a build directory:

..  code-block:: bash

    mkdir ${REPOSITORY_ROOT}/cmake_dir/ # create a subdirectory in the root of the repository
    cd ${REPOSITORY_ROOT}/cmake_dir/
    cmake .. # initialize the project

Cmake checks the dependencies and generates the main makefile that has three important targets:

    1. *molpher* -- This is the main target. It builds the C++ source code of the library and creates all necessary binary files. All targets mentioned below depend on this target.

    2. *molpher_install* -- This target will install the library in the given location.

        By default, this location is the :file:`dist/` folder in the repository root.
        This can be changed when the cmake project is initialized by setting
        `CMAKE_INSTALL_PREFIX <https://cmake.org/cmake/help/v3.5/variable/CMAKE_INSTALL_PREFIX.html>`_.
        Also by default, the TBB library is installed along with the molpher-lib files, because it
        is also a runtime dependency. If you do not want this behaviour you can instruct
        cmake to omit the installation of TBB with the :command:`-DINSTALL_TBB=OFF` option.

    3. *molpher_build_SWIG_Python* -- This will initialize the Python package along with the C++ extensions using the generated SWIG wrapper code.

        By default, the Python distribution under :command:`python` is used. You can specify a different executable
        by specifying the path with :command:`-DPYTHON_EXECUTABLE`
        when you run cmake.

        If you want to update the SWIG wrapping code before this target is run, you can instruct cmake to do so with
        the :command:`-DRUN_SWIG=ON` option.

        When this target finishes, all required files should be in place and you should be able to
        import the *molpher* Python package located under :file:`src/python/`. See the section below on how
        to install the package using *pip*.

        ..  warning:: As of yet, this target requires the *molpher_install* target above to install in the default location.

The build process of a given configuration is initialized from the cmake project directory as:

..  code-block:: bash

    make $CONFIG_NAME # CONFIG_NAME is the name of one of the configurations

Building the Python Extensions and Conda Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you built the library and initialized the Python package (the *molpher_build_SWIG_Python* target),
you can now simply install it with *pip*:

    ..  code-block:: bash

        pip install .

This will compile the Python extension and install the Python package to the currently active Python environment.

If you want to build the conda packages, you can a
python script located in the :file:`conda/` subdirectory of the repository root:

..  code-block:: bash

    python build_conda.py

..  attention:: You will need `conda-build <https://github.com/conda/conda-build>`_ and the  *jinja2* Python library to do that.

You can install the built packages by running:

..  code-block:: bash

    conda install --use-local molpher-lib

Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

This documentation is generated using the :file:`build_docs.sh`
script under the :file:`doc/` directory. However, in order to successfully build the documentation
you will need a few packages in your Python environment:

..  literalinclude:: ../../../doc/requirements.txt
    :language: none
    :caption: The requirements file for the environment used to build the documentation.

You can easily install all these packages with *pip*:

..  code-block:: bash

    cd ${REPOSITORY_ROOT}/doc/
    pip install -r requirements.txt

Also, make sure you initialized the Python package (built the *molpher_build_SWIG_Python* target as described above)
prior to generating the documentation so that Sphinx can successfully import the *molpher* package.

To update the GitHub pages with the current version of the documentation, you can simply
run the :file:`build_docs.sh` with the ```--upload``` option:

..  code-block:: bash

    build_docs.sh --upload

..  note:: You will need write access to the repository and an SSH key attached to your account to be able to upload
        the documentation.
