# Molpher-lib Changelog

Previous version: 0.0.0b0

## Changes
- made *cmake* the default build tool
- made the C++ code less verbose
- added pathfinder with antidecoys to examples

## Fixes
- prevented the extend operation from incrementing the distance
improvement counter of the newly appended leaves
- fixed a problem with thread and generation counts being reset when ExplorationTree.params was used to update morphing parameters
