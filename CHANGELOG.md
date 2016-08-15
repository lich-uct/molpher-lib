# Molpher-lib Changelog

Previous version: 0.0.0b0

## Changes
- made *cmake* the default build tool
- made the C++ code less verbose
- implemented a pathfinder with antidecoys
- added a script for building dependencies
- created the `molpher.algorithms` package to house implementations 
of complete pathfinding algorithms and moved the algorithms from the tutorial to it 
(it now includes the classic, bidirectional and antidecoys algorithm)
- it is now possible to change how morphs are sorted with a callback (`SortMorphsCallback`)
- new operation was added that will remove morphs from the list of candidates according to the information in candidates mask (`CleanMorphsOper`)
- created a separate cmake project for tests

## Fixes
- prevented the extend operation from incrementing the distance
improvement counter of the newly appended leaves
- fixed a problem with thread and generation counts being reset when `ExplorationTree.params` was used to update morphing parameters
- updated Boost to 1.50.0 to fix a bug that prevented a successful build of the *threads* library
