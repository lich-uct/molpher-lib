# Molpher-lib

This C++/Python library is based on the program [Molpher](https://github.com/siret/molpher) and implements systematic [chemical space](https://en.wikipedia.org/wiki/Chemical_space) exploration via a method called [molecular morphing](http://www.ncbi.nlm.nih.gov/pubmed/24655571). The method uses stochastic optimization to search for a 'path' in chemical space that connects two molecules by the means of small structural changes. It can be used to sample areas of chemical space that are biologically relevant and create a diverse set of novel potentially active compounds that can be analyzed further or used as a basis for conventional [virtual screening](https://en.wikipedia.org/wiki/Virtual_screening).

The purpose of the library is to make molecular morphing easily accessible to a wider range of developers. Because the current architecture of Molpher is fairly rigid and, in some cases, rather tightly coupled, it is not very developer-friendly. Furthermore, it requires anyone who wishes to even slightly modify the algorithm to delve deeply into its internals and recompile the whole backend code, which has a lot of dependencies. 

 The main focus during the development is on:

1. reducing the codebase to only the most essential parts and removing unnecessary dependencies and functionality that is not important for developers (only the most essential backend code was kept),
2. refactoring the code so that as many features as possible can be easily exposed via an API,
3. exposing the API of the library to a high-level programming language (such as Python) for ease of use.

An extensive documentation with a tutorial for all versions of the library is available [here](https://lich-uct.github.io/molpher-lib/).

## Installation

### Installation with Anaconda

At the moment the library is only built for the 64-bit Linux platform and distributed in the form of [conda packages hosted on Anaconda Cloud](https://anaconda.org/lich/molpher-lib). This is probably the easiest and preferred way to install the library and the associated Python package.

All you need to do is just either get the full [Anaconda](https://www.continuum.io/downloads) distribution or its lightweight variant, [Miniconda](http://conda.pydata.org/miniconda.html). It is essentially a Python distribution, package manager and virtual environment in one and makes setting up a development environment for your project very easy. After installing Anaconda/Miniconda all you need to do is to run the following in the Linux terminal:

```bash
conda install -c lich molpher-lib
```

This will automatically download the latest version of the library and the corresponding Python package and install everything in the currently active environment (for more information on environments and the `conda` command see [Conda Test Drive](http://conda.pydata.org/docs/test-drive.html).

And that's it! Now you can take a look at the [documentation](https://lich-uct.github.io/molpher-lib/) and do some molecular morphing.
