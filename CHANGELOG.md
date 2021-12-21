# Molpher-lib Changelog

Current Version: 0.0.0b4
Previous Version: 0.0.0b3

## New Features

- expose the SAScore threshold as a parameter of the exploration tree (`tree.params.sascoreMax`).
- implement a new chemical operator that allows specification of a chemical reaction with SMARTS

## Changes

- removed the original chemical operator enums -> operators will be set dynamically with pure instances in the future

## Fixes

-  use concurrent vector for candidates mask to improve thread safety
