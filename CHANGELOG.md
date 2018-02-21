# Molpher-lib Changelog

Previous version: 0.0.0b1

## Changes
- implemented a basic pathfinder for quick customized searches (`molpher.algorithms.pathfinders.BasicPathfinder`)
- some dependencies were updated (RDKit to v2017.09.1, Boost to v1.63)
- all dependencies are now linked dynamically and packaged separately during a conda build
- the tests are now included in the main CMakeLists file
- the MorphingOperator interface was implemented 
    - AddAtom operator class was added that adheres to this interface
    - RemoveAtom operator class was added that adheres to this interface
    - AddBond operator class was added that adheres to this interface
    - RemoveBond operator class was added that adheres to this interface
    - MutateAtom operator class was added that adheres to this interface
    - InterlayAtom operator class was added that adheres to this interface
    - ContractBond (originally BondContraction) operator class was added that adheres to this interface
    - RerouteBond (originally BondReroute) operator class was added that adheres to this interface
- MolpherMol can now be used to get all neighbors of an atom with the getNeighbors(int idx) method

## Fixes
- number of linked dependencies was reduced (unused libraries are not linked now)
