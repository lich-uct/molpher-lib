# Molpher-lib Changelog

Previous version: 0.0.0b2

## New Features

- the tree now has a `source` and `target` attributes which contain SMILES strings of the source and target structures

## Changes

- dependencies were updated to more recent versions, molpher-lib should now work fine with:
    - TBB 2019.9
    - BOOST 1.73.0
    - RDKIT 2020.09.1.0
- because RDKit is dropping Python 2 support (https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg08354.html), Python 2 conda packages will no longer be built for molpher-lib

## Fixes
- fixed issue #6 (Cryptic error raised when `max_iters` is specified in algorithm settings)
- fixes issue #9 (the newest versions of dependencies should work fine now)
