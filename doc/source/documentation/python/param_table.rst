..  csv-table:: Morphing parameters recognized by the current version.
    :header: "Attribute", "Default Value", "Brief Description", "Setter", "Getter"
    :name: param-table

    `source`, `None`, "SMILES of the :term:`source molecule`.", `setSource`, `getSource`
    `target`, `None`, "SMILES of the :term:`target molecule`.", `setTarget`, `getTarget`
    `operators`, "a `tuple` of :term:`selectors` [1]_", "A `tuple` of identifiers of the permitted `chemical operators`.", `setChemicalOperators`, `getChemicalOperators`
    `accept_max`, 100, "Maximum number of candidates accepted at once (based on their position in `ExplorationTree.candidates`).", `setCntCandidatesToKeepMax`, `getCntCandidatesToKeepMax`
    `accept_min`, 50, "Minimum number of candidates accepted during probability filtering.", `setCntCandidatesToKeep`, `getCntCandidatesToKeep`
    `close_produce`, 150, "Maximum number of morphs to produce with an `ExplorationTree.generateMorphs()` call when close to the :term:`target molecule`.", `setCntMorphsInDepth`, `getCntMorphsInDepth`
    `far_produce`, 80, "Maximum number of morphs to produce with an `ExplorationTree.generateMorphs()` call.", `setCntMorphs`, `getCntMorphs`
    `far_close_threshold`, 0.15, "Molecular distance below which the :term:`target molecule` and a :term:`morph` are considered to be close.", `setDistToTargetDepthSwitch`, `getDistToTargetDepthSwitch`
    `fingerprint`, `FP_MORGAN`, "Identification string of the current fingerprint strategy.", `setFingerprint`, `getFingerprint`
    `similarity`, "TANIMOTO", "Identification string of the current fingerprint strategy.", `setSimilarityCoefficient`, `getSimilarityCoefficient`
    `max_morphs_total`, 1500, "Maximum number of :term:`morphs <morph>` allowed to be derived from one molecule and the allowed number of non-producing descendants before a molecule is removed from the tree.", `setCntMaxMorphs`, `getCntMaxMorphs`
    `non_producing_survive`, 5, "Number of iterations before descendants of a non-producing molecule are removed from the tree.", `setItThreshold`, `getItThreshold`
    `weight_max`, 100000.0, "Maximum molecular weight of one morph.", `setMaxAcceptableMolecularWeight`, `getMaxAcceptableMolecularWeight`
    `weight_min`, 0.0, "Minimum molecular weight of one morph.", `setMinAcceptableMolecularWeight`, `getMinAcceptableMolecularWeight`

..  [1] ('ADD_ATOM', 'ADD_BOND', 'BOND_CONTRACTION', 'BOND_REROUTE', 'INTERLAY_ATOM', 'MUTATE_ATOM', 'REMOVE_ATOM', 'REMOVE_BOND')