
#include "molpher_API/MolpherMol.hpp"

MolpherMol::MolpherMol() : mol(new MolpherMolecule()), selfAllocated(true) {
    // no action
}


MolpherMol::MolpherMol(const MolpherMol& other) : mol(nullptr), selfAllocated(true) {
    MolpherMolecule* new_mol = new MolpherMolecule();
    *new_mol = *(other.mol);
    this->mol = new_mol;
}

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

MolpherMolecule& MolpherMol::fetchMolpherMolecule() const {
    return *mol;
}

bool MolpherMol::isBound() const {
    return !selfAllocated;
}

MolpherMol* MolpherMol::copy() const {
    return new MolpherMol(*this);
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

