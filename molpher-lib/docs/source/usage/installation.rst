Installing the Library
======================

Given the early stage of development, installing the library is a little bit tricky at the moment.
So far the library was only built and tested on a 64-bit Linux computer running the Debian
operating system. If you want to use the library with this (or similar) setup, you can just download
the current version (|version|) with all dependencies built from `here
<https://drive.google.com/file/d/0B2rizkCQQcoyaHJjaFFSdE9DdUk/view?usp=sharing>`_.

If you download this precompiled version, all you need to do is just unpack the linked archive
and you can use the library right away. If you need to install the Python package,
just do:

..  code-block:: bash

    cd molpher-lib # from the installation directory root
    ./build_python.sh
    pip install .

This will automatically build Python bindings and installs the library as a Python package into
the presently active Python environment.

Building the most recent version from source
--------------------------------------------

If you are a developer and want to build the most recent version from source,
you will need do everything from scratch. This section
describes how to do that on the Linux platform such as the one mentioned above
and hints on some possible problems one might encounter. If you want to build the library on
a different platform, this section might provide useful hints as well.

Prerequisites
~~~~~~~~~~~~~

Before you start building, there are a few prerequisites that need to be satisified on your system:

    - *build-essential* -- You will need this package to be able to build the software. It contains a compiler and other tools important for the build.
    - *g++-4.6* -- You will need this older version of g++ on your system as well. It is needed to build an old version of the RCF library. On Debian this package was moved to *oldstable*. Therefore, don't forget to add it to your :file:`/etc/apt/sources.list`. For example:

        ..  code-block:: bash

            deb http://ftp.cvut.cz/debian/ oldstable main non-free contrib
            deb-src http://ftp.cvut.cz/debian/ oldstable main non-free contrib

        ..  note:: There is effort to `remove or update
                in the future <https://github.com/siret/Molpher/tree/library_update>`_,
                but for the time being it is necessary to build it and it only can be done with older
                versions of g++.

    - *cmake* -- Used to build some of the dependencies.
    - *bison* and *flex* -- Required by the rdkit toolkit.
    - *setuptools* and *pip* -- If you want to use *pip* to install the Python package, you will need this to be installed with your Python distribution.
    - *python{version}-dev* -- You will need this package in order to build the Python bindings for the particular Python *version*.

When those requirements are satisfied you just have to checkout the latest version from Bitbucket:

..  code-block:: bash

    git clone -b molpher-lib https://{username}@bitbucket.org/Krwemil/molpher.git # (replace {username} with your Bitbucket username)

and start building.

..  note:: The repository is currently private and you will need the authorization of the author to access it.

Building Molpher and Its Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, download the `Molpher dependency bundle
<https://drive.google.com/file/d/0B2rizkCQQcoybFdhOFExaVk5c0U/view?usp=sharing>`_
and extract the contents to the :file:`dependecies/` directory of the repository root.

When it's done, execute the :file:`install-redist-linux-ia64.sh` script from the repository root.
This should automatically build all dependencies and compile Molpher as a static library
in the current directory hierarchy. If you want to change that, you can set the :makevar:`PREFIX`
variable in :file:`install-redist-linux-ia64.sh`. That will install everything to the specified directory.

Building the Molpher Library and Python Bindings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the `Molpher-backend project` and its dependencies have been built, we just need to compile
the `Molpher-lib project`. If you only need to use the C++ interface,
you can just do:

..  code-block:: bash

    cd molpher-lib
    make CONF=Linux64_Debug

If you want to build the Python bindings as well,
you can simply run:

..  code-block:: bash

    cd molpher-lib
    ./build_python.sh

To install as a package to the current Python do:

..  code-block:: bash

    cd molpher-lib
    pip install .

There is also a convenience script (``install_python.sh``) in the root of the repository
which does all of the above automatically.