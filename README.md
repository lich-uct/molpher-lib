# Molpher-lib: Introduction

This C++/Python library is a [chemical space](https://en.wikipedia.org/wiki/Chemical_space) exploration software. It is based on the [Molpher](https://github.com/siret/molpher) program which introduced a chemical space exploration method called [molecular morphing](http://www.ncbi.nlm.nih.gov/pubmed/24655571). The original Molpher method uses stochastic optimization to traverse chemical space between two existing molecules. The main promise of this algorithm is that a virtual library enriched in compounds with improved biological activity could be generated in this way. The library is based on the popular [RDKit](http://www.rdkit.org/) cheminformatics framework, which makes work with the generated structures very easy.

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
conda install -c conda-forge -c lich molpher-lib
```

This will automatically download the latest version of the library and install everything to the currently active environment (for more information on environments and the `conda` command see [Conda Test Drive](http://conda.pydata.org/docs/test-drive.html)). The library depends on a few other packages so do not forget to add the `conda-forge` channel.

If you are interested in the development snapshots of the library 
(most up to date code, but can contain bugs)
, use the `dev` channel instead:

```bash
conda install -c conda-forge -c lich/label/dev molpher-lib
```

After that the library should be importable in Python and you should be able to successfully run the integrated unit tests:

```python
from molpher.tests import run

run()
```

For a quick start, you can check out the [Jupyter notebooks](https://github.com/lich-uct/molpher-lib/tree/master/doc/notebooks) 
that contain examples from the [documentation](https://lich-uct.github.io/molpher-lib/latest/).

### Compiling from Source

Compiling and installing from source is a little bit more elaborate. This process is described in detail in the [documentation](https://lich-uct.github.io/molpher-lib/latest/usage/installation.html#building-and-installing-from-source-linux), 
but in the simplest case you can use the attached conda build environment 
to get all the dependencies and tools you will need for compilation:

```bash
# clone the repository
git clone https://github.com/lich-uct/molpher-lib.git
export REPOSITORY_ROOT=`pwd`/molpher-lib

# create the build environment from the attached file and activate it
cd ${REPOSITORY_ROOT}
git checkout dev # or the branch/tag/commit you want
conda env create -n {name} -f environment.yml # replace {name} with the name of your environment, i.e. molpher-lib-build
conda activate {name}

# build the library itself
mkdir cmake-build
cd cmake-build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DPYTHON_EXECUTABLE=python3 -DRUN_SWIG=ON # you can specify other options (see the documentation link above for more info)
make molpher_install_python # compile the library and Python wrappers, will also install to ${REPOSITORY_ROOT}/dist
```

#### Testing the Compiled Code

With the build environment still active, set the appropriate variables:

```bash
cd ${REPOSITORY_ROOT}
export PYTHONPATH=${REPOSITORY_ROOT}/dist/lib/python3.9/site-packages/molpher-0.0.0b3-py3.9-linux-x86_64.egg/ # change according to your Python and Molpher-lib version
export LD_LIBRARY_PATH=${REPOSITORY_ROOT}/dist/lib/:${CONDA_PREFIX}/lib/
```

and you should be able to import the built library from Python. You can verify the installation by running unit tests:

```bash
python3
```

Inside the interpreter:

```python
from molpher.tests import run

run()
```

If you want to explore some example code from the documentations, there are
a few Jupyter notebooks located under `doc/notebooks`. Just install Jupyter to your
build environment and set the variables as above and you can simply run the notebookes like so:

```bash
cd ${REPOSITORY_ROOT}/doc/notebooks/
jupyter-notebook
```

This installation process has been tested on common Debian-based systems so experience on other Linux flavors may differ. If you run into problems, report them to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues) and hopefully someone will be able to help.
