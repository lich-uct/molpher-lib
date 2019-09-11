# Molpher-lib: Introduction

This C++/Python library is a [chemical space](https://en.wikipedia.org/wiki/Chemical_space) exploration software. It is based on the [Molpher](https://github.com/siret/molpher) program which introduced a chemical space exploration method called [molecular morphing](http://www.ncbi.nlm.nih.gov/pubmed/24655571). The original Molpher method uses stochastic optimization to traverse chemical space between two existing molecules. The main promise of this algorithm is that a virtual library enriched in compounds with improved biological activity could be generated in this way.

The purpose of Molpher-lib is to bring molecular morphing closer to the cheminformatics community, but also offer new features that go beyond the capabilities of the original Molpher program. Molpher-lib makes it possible to roam the chemical universe freely and with little constraints on the inputs. For example, we could just use a carbon atom as a starting point and have Molpher-lib autonomously evolve it into a complete molecular structure. To ensure that the generated molecules have required properties, Molpher-lib also helps with implementation of custom rules and constraints. If you want to know more about Molpher-lib and its usage, make sure to check out some [examples on the website](https://lich-uct.github.io/molpher-lib/examples.html). We also have some [Jupyter notebooks](https://github.com/lich-uct/molpher-lib/tree/master/doc/notebooks) with examples 
that you can explore.

If you would like to participate in the development or just check out the current features of the library, there is extensive [documentation](https://lich-uct.github.io/molpher-lib/latest/) which can help you. A big part of the documentation is dedicated to a detailed [tutorial](https://lich-uct.github.io/molpher-lib/latest/usage/tutorial.html) that should introduce the philosophy of Molpher-lib in more detail and give you a good idea of what it is currently capable of. 

The library is actively developed and many new features are planned to be added. The long-term goal is to make Molpher-lib a universal and easy-to-use *de novo* drug design framework. Ideas, comments and feature requests are more than welcome and can be submitted to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues). You can also [subscribe](https://github.com/lich-uct/molpher-lib/commits/dev.atom) to the RSS feed of the [dev](https://github.com/lich-uct/molpher-lib/tree/dev) branch for development updates. If you want to know what is new in the current version, you can look at the [changelog](CHANGELOG.md).

## Installation

Supported platforms:

  - Linux 64-bit

At the moment, the library binaries are only compiled for 64-bit Linux systems. However, development for other platforms is also planned. If you manage to compile the library on a different platform, consider making a pull request or comment on the [issue tracker](https://github.com/lich-uct/molpher-lib/issues). Any help is much appreciated.

### Installation with Anaconda

Molpher-lib is distributed as a [conda package](https://anaconda.org/lich/molpher-lib). At the moment, this is the preferred way to install and use the library. All you need to do is get the full [Anaconda](https://www.continuum.io/downloads) distribution or its lightweight variant, [Miniconda](http://conda.pydata.org/miniconda.html). It is essentially a Python distribution, package manager and virtual environment in one and makes setting up a development environment for any project very easy. After installing Anaconda/Miniconda you can run the following in the Linux terminal:

```bash
conda install -c rdkit -c lich molpher-lib
```

This will automatically download the latest version of the library and install everything to the currently active environment (for more information on environments and the `conda` command see [Conda Test Drive](http://conda.pydata.org/docs/test-drive.html)). The library depends on the popular cheminformatics toolkit [RDKit](http://rdkit.org) so do not forget to add the rdkit channel.

If you are interested in the development snapshots of the library 
(most up to date code, but can contain bugs)
, you can use the `dev` channel instead:

```bash
conda install -c rdkit -c lich/label/dev molpher-lib
```

After that the library should import in your environment and you should be able to successfully run the integrated unit tests:

```python
from molpher.tests import run

run()
```

You can also check out the [Jupyter notebooks](https://github.com/lich-uct/molpher-lib/tree/master/doc/notebooks) with examples from the [documentation](https://lich-uct.github.io/molpher-lib/latest/).

### Compiling from Source

Compiling and installing from source is a little bit more elaborate. This process is described in detail in the [documentation](https://lich-uct.github.io/molpher-lib/latest/usage/installation.html#building-and-installing-from-source-linux), but in the simplest case the following should work:

```bash
# get dependencies
sudo apt-get install git build-essential python3-dev python3-numpy cmake python3-setuptools

# clone the repo
git clone https://github.com/lich-uct/molpher-lib.git
git checkout dev # or the branch/tag/commit you want
REPOSITORY_ROOT=`pwd`/molpher-lib

# this might take a while, but you if you are lucky, 
# cmake might be able to find dependencies 
# if you already have them somewhere on your system
# so you can skip this step if you have TBB, Boost and RDKit
# installed at standard locations
cd ${REPOSITORY_ROOT}/deps
./build_deps.sh --all

# finally, build the library itself
cd ${REPOSITORY_ROOT}
mkdir cmake-build
cd cmake-build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DPYTHON_EXECUTABLE=python3
make molpher_install_python
```

After setting the appropriate variables:

```bash
export CMAKE_INSTALL_PREFIX="${REPOSITORY_ROOT}/dist"
export DEPS_DIR=${CMAKE_INSTALL_PREFIX}/../deps
export PYTHONPATH=${DEPS_DIR}/rdkit/:${CMAKE_INSTALL_PREFIX}/lib/python3.5/site-packages
export LD_LIBRARY_PATH=${DEPS_DIR}/tbb/lib/intel64/gcc4.7:${DEPS_DIR}/rdkit/lib/:${DEPS_DIR}/boost/stage/lib:${CMAKE_INSTALL_PREFIX}/lib
```

you should be good to go:

```bash
python3
```

```python
from molpher.tests import run

run()
```

This will run the integrated unit tests. They should all pass without problems.

If you want to explore some example code from the documentations, there are
a few Jupyter notebooks located under `doc/notebooks`. You can create 
the needed conda environment 
and launch your Jupyter server as follows:

```bash
cd ${REPOSITORY_ROOT}
conda env create -f "environment.yml"
. source_2_activate
python setup.py build_ext --inplace
cd doc/notebooks/
jupyter-notebook
```

Note that you will need to have the library already compiled and installed in the standard 
`${REPOSITORY_ROOT}/dist` directory.

This installation process has been tested on common Debian-based systems so experience on other Linux flavors may differ. If you run into problems, report them to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues) and hopefully someone will be able to help.
