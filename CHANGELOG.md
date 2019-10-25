# Molpher-lib Changelog

Previous version: 0.0.0b2

## New Features

- the tree now has a `source` and `target` attributes which contain SMILES strings of the source and target structures

## Changes

- dependencies were updated to more recent versions, molpher-lib should now work fine with:
    - TBB 2019_U8
    - BOOST 1.67.0
    - RDKIT 2019.03.4
    - newer or older versions of the packages above could also work, but this was not tested
- because RDKit is dropping Python 2 support (https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg08354.html), Python 2 conda packages will no longer be built for molpher-lib as well and Python 2 will also not be supported in the project anymore
- `asRDMol` method of `MolpherMol` in C++ now returns the actual pointer of the underlying instance and not a copy anymore

## Fixes
- fixed issue #6 (Cryptic error raised when `max_iters` is specified in algorithm settings)
