# Molpher-lib Changelog

Previous version: 0.0.0b1

## Changes
- implemented a basic pathfinder for quick customized searches (`molpher.algorithms.pathfinders.BasicPathfinder`)
- dependencies were updated (the code should work fine with RDKit 2018.03.1, Boost 1.65 and TBB 2018 U3)
- all dependencies are now linked dynamically 
- conda environment libraries can now be searched during cmake configuration by setting the `CONDA_PREFIX` environment variable
- conda packages now depend on tbb, rdkit and boost libraries from anaconda repositories both during build time and runtime
- the unit tests binary is now included in the `CMakeLists.txt` file and it can also be installed with the project
- SWIG was updated to version 3.0.12
- other small optimizations of the build process
- the `MorphingOperator` interface was implemented and tested in both Python and C++
    - `AddAtom` operator class was added that adheres to this interface
    - `RemoveAtom` operator class was added that adheres to this interface
    - `AddBond` operator class was added that adheres to this interface
    - `RemoveBond` operator class was added that adheres to this interface
    - `MutateAtom` operator class was added that adheres to this interface
    - `InterlayAtom` operator class was added that adheres to this interface
    - `ContractBond` (originally BondContraction) operator class was added that adheres to this interface
    - `RerouteBond` (originally BondReroute) operator class was added that adheres to this interface
- the `Molpher` class was created and can now be used to aggregate operators and generate random morphs of a structure (can use multiple threads)
- `ExplorationTree` and `GenerateMorphsOper` operation now use the new `Molpher` interface to facilitate exploration
- `ExplorationTree` and `GenerateMorphsOper` are also able to take instances of `MorphCollector`, which can be used to explore and/or modify morphs right after they are generated
- `MolpherMol` can be used to get all neighbors of an atom with the `getNeighbors(int idx)` method
- `MolpherMol` can be initialized from an MDL MOL block (using the `fromMolBlock()` static method) and from and RDKit molecule (using the `other` parameter in its constructor)
- `MolpherMol` can readily be converted to an RDKit molecule with the `asRDMol()` method
- `ExplorationTree` has a new `fetchPathTo` method, which backtracks from the given structure and returns the sequence of `MolpherMol` instances that lead to its creation, starting from the source molecule.

## Fixes
- number of linked dependencies was reduced and unused libraries are now not linked
