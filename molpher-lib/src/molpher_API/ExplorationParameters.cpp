
#include "molpher_API/ExplorationParameters.hpp"

ExplorationParameters::ExplorationParameters() {
    // no action
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

MolpherMol ExplorationParameters::getSourceMol() const {
    return MolpherMol(iterSnapshot.source.smile);
}
