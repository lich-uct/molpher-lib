# Molpher-lib Changelog

Current Version: 0.0.0b3
Previous Version: 0.0.0b2

## New Features

- the tree now has a `source` and `target` attributes which contain SMILES strings of the source and target structures

## Changes

- the `./deps/build_deps.sh` script is now retired and removed
  - use the conda build environment instead (its current snapshot is always saved in `environment.yml`)
  - the environment can be supplied to cmake via the CONDA_PREFIX environment variable (automatically set when `conda activate {env}` is run)
- because [RDKit dropped Python 2 support](https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg08354.html), Python 2 conda packages will no longer be built for molpher-lib
- simplified conda package build
- use `conda-forge` as the main channel for dependencies (also necessary to specify during `conda install`), README changed accordingly
- the atom types incorporated into generated structures now follow realistic distribution derived from ChEMBL
  - this affect the `AddAtom`, `InterlayAtom` and `MutateAtom` operators 

## Fixes
- fixed issue #6 (Cryptic error raised when `max_iters` is specified in algorithm settings)
- fixes issue #9 (the newest versions of dependencies should work fine now)
- fixed a few minor issues that caused some dependencies to break when certain dependency versions were used
