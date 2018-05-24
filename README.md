# Molpher-lib

This C++/Python library is a [chemical space](https://en.wikipedia.org/wiki/Chemical_space) exploration software. It is based on the [Molpher](https://github.com/siret/molpher) program which implements a method called [molecular morphing](http://www.ncbi.nlm.nih.gov/pubmed/24655571). This method uses stochastic optimization to traverse chemical space between two molecules. It can be used to sample unexplored areas that might contain new bioactive compounds with increased probability. The purpose of the library is to make molecular morphing more accessible and flexible and to provide good basis for further experimentation in this area. See the [official website](https://lich-uct.github.io/molpher-lib/) for additional information and usage examples.

The library is actively developed and a lot of new features are planned for the future. The long-term goal is to make Molpher-lib a universal and easy-to-use *de novo* drug design framework with possibilities that go beyond molecular morphing. If this seems interesting to you, you can take a look at the [documentation](https://lich-uct.github.io/molpher-lib/latest/) to get an idea of what the library is currently capable of as well as what features are planned for future releases. Ideas, comments or feature requests are more than welcome and can be submitted to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues). You can also [subscribe](https://github.com/lich-uct/molpher-lib/commits/master.atom) to the RSS feed for updates. If you want to know what is new in the current version, you can look at the [changelog](CHANGELOG.md).

At the moment, the library is only intended for use on 64-bit Linux systems. However, development for other platforms is also a priority. If you manage to compile the library on a different platform, consider making a pull request or comment on the [issue tracker](https://github.com/lich-uct/molpher-lib/issues). Any help is much appreciated.

## Installation

### Installation with Anaconda

Molpher-lib is distributed as a [conda package](https://anaconda.org/lich/molpher-lib). At the moment, this is the preferred way to install and use the library. All you need to do is just either get the full [Anaconda](https://www.continuum.io/downloads) distribution or its lightweight variant, [Miniconda](http://conda.pydata.org/miniconda.html). It is essentially a Python distribution, package manager and virtual environment in one and makes setting up a development environment for your project very easy. After installing Anaconda/Miniconda you can run the following in the Linux terminal:

```bash
conda install -c rdkit -c lich molpher-lib
```

This will automatically download the latest version of the library and install everything to the currently active environment (for more information on environments and the `conda` command see [Conda Test Drive](http://conda.pydata.org/docs/test-drive.html)).

### Installation from Source

Installing from source is a little bit more elaborate because Molpher-lib contains a lot of C++ code that needs to be compiled first. This process is described [here](https://lich-uct.github.io/molpher-lib/latest/usage/installation.html#building-and-installing-from-source-linux) in detail, but in the simplest case the following should work:

```bash
git clone https://github.com/lich-uct/molpher-lib.git
REPOSITORY_ROOT=`pwd`/molpher-lib

# this might take a while, but you if you are lucky, 
# cmake might be able to find dependencies 
# if you already have them somewhere on your system
# so you can skip this step if you have TBB, Boost and RDKit
# installed at standard locations on your system
cd ${REPOSITORY_ROOT}/deps
./build_deps.sh --all

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

```python
from molpher.tests import run

run()
```

This installation process was only tested on common Debian-based systems so experience on other Linux flavors may be different. If you run into problems, report them to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues) and hopefully someone will be able to help.
