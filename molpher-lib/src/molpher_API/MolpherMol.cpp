
#include "molpher_API/MolpherMol.hpp"

MolpherMol::MolpherMol() : mol(new MolpherMolecule()), selfAllocated(true) {
    // no action
}


MolpherMol::MolpherMol(const std::string& smile) : selfAllocated(true) {
    std::string temp(smile);
    mol = new MolpherMolecule(temp);
}

MolpherMol::MolpherMol(MolpherMolecule& mol) : selfAllocated(false) {
    this->mol = &mol;
}

MolpherMol::~MolpherMol() {
    if (selfAllocated) delete mol;
}

MolpherMolecule& MolpherMol::getMol() {
    return *mol;
}

std::string MolpherMol::getSMILES() {
    return mol->smile;
}

double MolpherMol::getDistToTarget() {
    return mol->distToTarget;
}

std::string MolpherMol::getParentSMILES() {
    return mol->parentSmile;
}

