..  csv-table:: Morphing parameters recognized by the current implementation.
    :header: "Attribute", "Setter", "Getter", "Default Value", "Brief Description"
    :name: param-table

    `source`, `setSource`, `getSource`, "''", "SMILES of the `source molecule`."
    `target`, `setTarget`, `getTarget`, "''", "SMILES of the `target molecule`."
    `operators`, `setChemicalOperators`, `getChemicalOperators`, "('ADD_ATOM', 'ADD_BOND', 'BOND_CONTRACTION', 'BOND_REROUTE', 'INTERLAY_ATOM', 'MUTATE_ATOM', 'REMOVE_ATOM', 'REMOVE_BOND')", "A `tuple` of identifiers of the permitted `chemical operators`."
    `accept_max`, `setCntCandidatesToKeepMax`, `getCntCandidatesToKeepMax`, 100, "Maximum number of candidates accepted at once (based on their position in `ExplorationTree.candidates`)."
    `accept_min`, `setCntCandidatesToKeep`, `getCntCandidatesToKeep`, 50, "Minimum number of candidates accepted during probability filtering."
    `close_produce`, `setCntMorphsInDepth`, `getCntMorphsInDepth`, 150, "Maximum number of morphs to produce with an `ExplorationTree.generateMorphs()` call when close to the `target molecule`."
    `far_produce`, `setCntMorphs`, `getCntMorphs`, 80, "Maximum number of morphs to produce with an `ExplorationTree.generateMorphs()` call."
    `far_close_threshold`, `setDistToTargetDepthSwitch`, `getDistToTargetDepthSwitch`, 0.15, "Molecular distance below which the `target molecule` and a `morph` are cosidered to be close."
    `fingerprint`, `setFingerprint`, `getFingerprint`, "MORGAN", "Identification string of the current fingerprint strategy."
    `similarity`, `setSimilarityCoefficient`, `getSimilarityCoefficient`, "TANIMOTO", "Identification string of the current fingerprint strategy."
    `max_morphs_total`, `setCntMaxMorphs`, `getCntMaxMorphs`, 1500, "Maximum number of morphs allowed to be derived from one molecule and the allowed number of non-producing descendents before a molecule is removed from the tree."
    `non_producing_survive`, `setItThreshold`, `getItThreshold`, 5, "Number of iterations before descendents of a non-producing molecule are removed from the tree."
    `weight_max`, `setMaxAcceptableMolecularWeight`, `getMaxAcceptableMolecularWeight`, 100000.0, "Maximum molecular weight of one morph."
    `weight_min`, `setMinAcceptableMolecularWeight`, `getMinAcceptableMolecularWeight`, 0.0, "Minimum molecular weight of one morph."