
#include <algorithm>

#include "ExplorationParameters.hpp"

std::map<std::string, ChemOperSelector> ExplorationParameters::CHEMOPERS_SELECTOR_MAPPPING {
        std::make_pair("ADD_ATOM", ChemOperSelector::OP_ADD_ATOM),
        std::make_pair("REMOVE_ATOM", ChemOperSelector::OP_REMOVE_ATOM),
        std::make_pair("ADD_BOND", ChemOperSelector::OP_ADD_BOND),
        std::make_pair("REMOVE_BOND", ChemOperSelector::OP_REMOVE_BOND),
        std::make_pair("MUTATE_ATOM", ChemOperSelector::OP_MUTATE_ATOM),
        std::make_pair("INTERLAY_ATOM", ChemOperSelector::OP_INTERLAY_ATOM),
        std::make_pair("BOND_REROUTE", ChemOperSelector::OP_BOND_REROUTE),
        std::make_pair("BOND_CONTRACTION", ChemOperSelector::OP_BOND_CONTRACTION)
    };

std::map<std::string, SimCoeffSelector> ExplorationParameters::SIMILARITY_SELECTOR_MAPPPING {
    std::make_pair("ALL_BIT", SimCoeffSelector::SC_ALL_BIT),
    std::make_pair("ASYMMETRIC", SimCoeffSelector::SC_ASYMMETRIC),
    std::make_pair("BRAUN_BLANQUET", SimCoeffSelector::SC_BRAUN_BLANQUET),
    std::make_pair("COSINE", SimCoeffSelector::SC_COSINE),
    std::make_pair("DICE", SimCoeffSelector::SC_DICE),
    std::make_pair("KULCZYNSKI", SimCoeffSelector::SC_KULCZYNSKI),
    std::make_pair("MC_CONNAUGHEY", SimCoeffSelector::SC_MC_CONNAUGHEY),
    std::make_pair("ON_BIT", SimCoeffSelector::SC_ON_BIT),
    std::make_pair("RUSSEL", SimCoeffSelector::SC_RUSSEL),
    std::make_pair("SOKAL", SimCoeffSelector::SC_SOKAL),
    std::make_pair("TANIMOTO", SimCoeffSelector::SC_TANIMOTO),
    std::make_pair("TVERSKY_SUBSTRUCTURE", SimCoeffSelector::SC_TVERSKY_SUBSTRUCTURE),
    std::make_pair("TVERSKY_SUPERSTRUCTURE", SimCoeffSelector::SC_TVERSKY_SUPERSTRUCTURE)
};

std::map<std::string, FingerprintSelector> ExplorationParameters::FINGERPRINT_SELECTOR_MAPPPING {
    std::make_pair("ATOM_PAIRS", FingerprintSelector::FP_ATOM_PAIRS),
    std::make_pair("MORGAN", FingerprintSelector::FP_MORGAN),
    std::make_pair("TOPOLOGICAL", FingerprintSelector::FP_TOPOLOGICAL),
    std::make_pair("TOPOLOGICAL_LAYERED_1", FingerprintSelector::FP_TOPOLOGICAL_LAYERED_1),
    std::make_pair("TOPOLOGICAL_LAYERED_2", FingerprintSelector::FP_TOPOLOGICAL_LAYERED_2),
    std::make_pair("VECTORFP", FingerprintSelector::FP_VECTORFP),
    std::make_pair("TOPOLOGICAL_TORSION", FingerprintSelector::FP_TOPOLOGICAL_TORSION),
    std::make_pair("EXT_ATOM_PAIRS", FingerprintSelector::FP_EXT_ATOM_PAIRS),
    std::make_pair("EXT_MORGAN", FingerprintSelector::FP_EXT_MORGAN),
    std::make_pair("EXT_TOPOLOGICAL", FingerprintSelector::FP_EXT_TOPOLOGICAL),
    std::make_pair("EXT_TOPOLOGICAL_LAYERED_1", FingerprintSelector::FP_EXT_TOPOLOGICAL_LAYERED_1),
    std::make_pair("EXT_TOPOLOGICAL_LAYERED_2", FingerprintSelector::FP_EXT_TOPOLOGICAL_LAYERED_2),
    std::make_pair("EXT_TOPOLOGICAL_TORSION", FingerprintSelector::FP_EXT_TOPOLOGICAL_TORSION),
};
    
ExplorationParameters::ExplorationParameters() {
    // no action
}

ExplorationParameters::ExplorationParameters(const ExplorationParameters& other) {
    iterSnapshot = other.iterSnapshot;
}

//IterationSnapshot ExplorationParameters::createIterationSnapshot() const {
//    return iterSnapshot;
//}

bool ExplorationParameters::valid() {
    bool validity = false;
    if (iterSnapshot.target.smile.empty()) {
        iterSnapshot.target.smile = "TARGET NOT SET";
        validity = iterSnapshot.IsValid();
        iterSnapshot.target.smile = "";
    } else {
        validity = iterSnapshot.IsValid();
    }
    return validity;
}

void ExplorationParameters::setSourceMol(const std::string& mol) {
    std::string temp(mol);
    iterSnapshot.source = MolpherMolecule(temp);
}

void ExplorationParameters::setSourceMol(MolpherMol& mol) {
    iterSnapshot.source = mol.fetchMolpherMolecule();
}

MolpherMol* ExplorationParameters::getSourceMol() {
    return new MolpherMol(iterSnapshot.source, true);
}

MolpherMol* ExplorationParameters::getTargetMol() {
    return new MolpherMol(iterSnapshot.target, true);
}

void ExplorationParameters::setTargetMol(const std::string& mol) {
    std::string temp(mol);
    iterSnapshot.target = MolpherMolecule(temp);
}

void ExplorationParameters::setTargetMol(MolpherMol& mol) {
    iterSnapshot.target = mol.fetchMolpherMolecule();
}

const std::vector<std::string>& ExplorationParameters::getChemOperators() {
    auto ret = new std::vector<std::string>();
    for (auto& pair : CHEMOPERS_SELECTOR_MAPPPING) {
        if (std::find(iterSnapshot.chemOperSelectors.begin(), iterSnapshot.chemOperSelectors.end(), pair.second) != iterSnapshot.chemOperSelectors.end()) {
            ret->push_back(pair.first);
        }
    }
    return *ret;
}

void ExplorationParameters::setChemOperators(const std::vector<std::string>& choices) {
    iterSnapshot.chemOperSelectors.clear();
    for (auto& choice : choices) {
        if (CHEMOPERS_SELECTOR_MAPPPING.find(choice) != CHEMOPERS_SELECTOR_MAPPPING.end()) {
            iterSnapshot.chemOperSelectors.push_back(CHEMOPERS_SELECTOR_MAPPPING[choice]);
        } else {
            throw std::runtime_error("Unknown morphing operator: " + choice);
        }
    }
}

int ExplorationParameters::getCntCandidatesToKeep() {
    return iterSnapshot.params.cntCandidatesToKeep;
}

int ExplorationParameters::getCntCandidatesToKeepMax() {
    return iterSnapshot.params.cntCandidatesToKeepMax;
}

int ExplorationParameters::getCntMaxMorphs() {
    return iterSnapshot.params.cntMaxMorphs;
}

int ExplorationParameters::getCntMorphs() {
    return iterSnapshot.params.cntMorphs;
}

int ExplorationParameters::getCntMorphsInDepth() {
    return iterSnapshot.params.cntMorphsInDepth;
}

double ExplorationParameters::getDistToTargetDepthSwitch() {
    return iterSnapshot.params.distToTargetDepthSwitch;
}

int ExplorationParameters::getItThreshold() {
    return iterSnapshot.params.itThreshold;
}

double ExplorationParameters::getMaxAcceptableMolecularWeight() {
    return iterSnapshot.params.maxAcceptableMolecularWeight;
}

double ExplorationParameters::getMinAcceptableMolecularWeight() {
    return iterSnapshot.params.minAcceptableMolecularWeight;
}

void ExplorationParameters::setCntCandidatesToKeep(int value) {
    iterSnapshot.params.cntCandidatesToKeep = value;
}

void ExplorationParameters::setCntCandidatesToKeepMax(int value) {
    iterSnapshot.params.cntCandidatesToKeepMax = value;
}

void ExplorationParameters::setCntMaxMorphs(int value) {
    iterSnapshot.params.cntMaxMorphs = value;
}

void ExplorationParameters::setCntMorphs(int value) {
    iterSnapshot.params.cntMorphs = value;
}

void ExplorationParameters::setCntMorphsInDepth(int value) {
    iterSnapshot.params.cntMorphsInDepth = value;
}

void ExplorationParameters::setDistToTargetDepthSwitch(double value) {
    iterSnapshot.params.distToTargetDepthSwitch = value;
}

void ExplorationParameters::setItThreshold(int value) {
    iterSnapshot.params.itThreshold = value;
}

void ExplorationParameters::setMaxAcceptableMolecularWeight(double weight) {
    iterSnapshot.params.maxAcceptableMolecularWeight = weight;
}

void ExplorationParameters::setMinAcceptableMolecularWeight(double weight) {
    iterSnapshot.params.minAcceptableMolecularWeight = weight;
}

std::string ExplorationParameters::getFingerprint() {
    for (auto& fp : FINGERPRINT_SELECTOR_MAPPPING) {
        if (fp.second == iterSnapshot.fingerprintSelector) {
            return fp.first;
        }
    }
    throw std::runtime_error("No matching fingerprint selector found.");
}

std::string ExplorationParameters::getSimilarityCoef() {
    for (auto& param : SIMILARITY_SELECTOR_MAPPPING) {
        if (param.second == iterSnapshot.simCoeffSelector) {
            return param.first;
        }
    }
    throw std::runtime_error("No matching similarity coefficient selector found.");
}

void ExplorationParameters::setFingerprint(const std::string& fp) {
    if (FINGERPRINT_SELECTOR_MAPPPING.find(fp) != FINGERPRINT_SELECTOR_MAPPPING.end()) {
        iterSnapshot.fingerprintSelector = FINGERPRINT_SELECTOR_MAPPPING[fp];
    } else {
        throw std::runtime_error("No matching fingerprint selector found for: " + fp);
    }
}

void ExplorationParameters::setSimilarityCoef(const std::string& coef) {
    if (SIMILARITY_SELECTOR_MAPPPING.find(coef) != SIMILARITY_SELECTOR_MAPPPING.end()) {
        iterSnapshot.simCoeffSelector = SIMILARITY_SELECTOR_MAPPPING[coef];
    } else {
        throw std::runtime_error("No matching similarity coefficient selector found for: " + coef);
    }
}



