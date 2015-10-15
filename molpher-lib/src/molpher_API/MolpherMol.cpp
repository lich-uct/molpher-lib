
#include "molpher_API/MolpherMol.hpp"


MolpherMol::MolpherMol(const std::string& smile) {
    std::string temp(smile);
    mol = MolpherMolecule(temp);
}

MolpherMol::MolpherMol(const MolpherMolecule& mol) {
    this->mol = mol;
}

MolpherMolecule& MolpherMol::getMol() {
    return mol;
}