
#include "molpher_API/MolpherMol.hpp"

MolpherMol::MolpherMol() : mol(new MolpherMolecule()), selfAllocated(true) {
    // no action
}


//MolpherMol::MolpherMol(const std::string& smile) : selfAllocated(true) {
//    std::string temp(smile);
//    mol = new MolpherMolecule(temp);
//}

MolpherMol::MolpherMol(MolpherMolecule& mol) : mol(&mol), selfAllocated(false) {
    // no action
}

MolpherMol::~MolpherMol() {
    if (selfAllocated) delete mol;
}

MolpherMol& MolpherMol::operator=(const MolpherMol& other) {
    MolpherMolecule* new_mol = new MolpherMolecule();
    *new_mol = *(other.mol);
    if (selfAllocated) {
        delete this->mol;
        this->mol = nullptr;
    }
    this->mol = new_mol;
    this->selfAllocated = true;
    return *this;
}

MolpherMolecule& MolpherMol::getMol() const {
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

