# Molpher-lib

This C++/Python library is a [chemical space](https://en.wikipedia.org/wiki/Chemical_space) exploration software. It is based on the [Molpher](https://github.com/siret/molpher) program which implements a method called [molecular morphing](http://www.ncbi.nlm.nih.gov/pubmed/24655571). This method uses stochastic optimization to traverse chemical space between two molecules. It can be used to sample unexplored areas that might contain new bioactive compounds with increased probability. The purpose of the library is to make molecular morphing more accessible and flexible and to provide good basis for further experimentation in this area. See the [official website](https://lich-uct.github.io/molpher-lib/) for additional information and usage examples.

The library is actively developed and a lot of new features are planned for the future. The long-term goal is to make Molpher-lib a universal and easy-to-use *de novo* drug design framework with possibilities that go beyond molecular morphing. If this seems interesting to you, you can take a look at the [documentation](https://lich-uct.github.io/molpher-lib/latest/) to get an idea of what the library is currently capable of as well as what features are planned for future releases. Ideas, comments or feature requests are more than welcome and can be submitted to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues). You can also [subscribe](https://github.com/lich-uct/molpher-lib/commits/master.atom) to the RSS feed for updates. If you want to know what is new in the current version, you can look at the [changelog](CHANGELOG.md).

At the moment, the library is only intended for use on 64-bit Linux systems. However, development for other platforms is also a priority. If you manage to compile the library on a different platform, consider making a pull request or comment on the [issue tracker](https://github.com/lich-uct/molpher-lib/issues). Any help is much appreciated.

## Installation

### Installation with Anaconda

Molpher-lib is distributed as a [conda package](https://anaconda.org/lich/molpher-lib). At the moment, this is the preferred way to install and use the library. All you need to do is just either get the full [Anaconda](https://www.continuum.io/downloads) distribution or its lightweight variant, [Miniconda](http://conda.pydata.org/miniconda.html). It is essentially a Python distribution, package manager and virtual environment in one and makes setting up a development environment for your project very easy. After installing Anaconda/Miniconda you can run the following in the Linux terminal:

```bash
conda install -c lich molpher-lib
```

This will automatically download the latest version of the library and install everything to the currently active environment (for more information on environments and the `conda` command see [Conda Test Drive](http://conda.pydata.org/docs/test-drive.html)).

### Installation from Source

Installing from source is a little bit more elaborate because Molpher-lib contains a lot of C++ code that needs to be compiled first. This process is described [here](https://lich-uct.github.io/molpher-lib/latest/usage/installation.html#building-and-installing-from-source-linux) in detail, but in the simplest case the following should work:

```bash
git clone https://github.com/lich-uct/molpher-lib.git
ROOT_DIR=`pwd`/molpher-lib
cd $ROOT_DIR/deps
./build_deps.sh --all # this might take a while, but you can bypass this if you already have Boost and RDKit compiled somewhere (see https://lich-uct.github.io/molpher-lib/)
cd $ROOT_DIR
mkdir build # the name of the directory does not matter here
cd build
cmake ..
make molpher_install_python

# optionally the python package can be tested:
cd $ROOT_DIR
python setup.py test
```

The make target above builds the library and installs everything  to `$ROOT_DIR/dist`. This is the default value for the `CMAKE_INSTALL_PREFIX` variable and it can be changed doing `cmake .. -DCMAKE_INSTALL_PREFIX=custom/install/directory/` instead of just plain `cmake ..`. This folder can be anywhere on the user's system provided that the following variables are set during runtime:

```bash
export PYTHONPATH=$CMAKE_INSTALL_PREFIX/lib/pythonX.Y/site-packages # replace X.Y with your Python version
export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib
```

The `molpher` package should now be importable from Python.

If you want to use the Python package right after the build, you can do so by just adding the `$ROOT_DIR/src/python` folder to PYTHONPATH like so:

```bash
export PYTHONPATH=$ROOT_DIR/src/python
```

This installation process was only tested on Debian 8.5 so experience on other Linux flavors may be different. If you run into problems, report them to the [issue tracker](https://github.com/lich-uct/molpher-lib/issues) and hopefully someone will be able to help.
