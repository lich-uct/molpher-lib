..  csv-table:: Filtering options and their descriptions.
    :header: "Option", "Description"

    `PROBABILITY`, "Filter the `candidate morphs` according to probability based on their position in `ExplorationTree.candidates`."
    `WEIGHT`, "Remove molecules that do not satisfy the weight constraints."
    `SYNTHESIS`, "The inbuilt synthetic feasibility filter."
    `MAX_DERIVATIONS`, "Maximum number of morphs derived from one molecule."
    `DUPLICATES`, "Remove candidates that are already in the tree (always on)."
    `HISTORIC_DESCENDENTS`, "Remove candidates that have already been tried by their parents."
    `ALL`, "Use all of the above filters."
