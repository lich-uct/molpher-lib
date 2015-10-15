
#include "molpher_API/ExplorationParameters.hpp"

ExplorationParameters::ExplorationParameters() {
    iterSnapshot.target.smile = "TARGET NOT SET";
}

//IterationSnapshot ExplorationParameters::createIterationSnapshot() const {
//    return iterSnapshot;
//}

bool ExplorationParameters::valid() const {
    return iterSnapshot.IsValid();
}

void ExplorationParameters::setSourceMol(const std::string& mol) {
    std::string temp(mol);
    iterSnapshot.source = MolpherMolecule(temp);
}

MolpherMol ExplorationParameters::getSourceMol() const {
    return MolpherMol(iterSnapshot.source.smile);
}
