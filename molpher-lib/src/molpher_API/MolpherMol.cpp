
#include "molpher_API/MolpherMol.hpp"


MolpherMol::MolpherMol(const std::string& smile) {
    std::string temp(smile);
    mol = MolpherMolecule(temp);
}

MolpherMolecule& MolpherMol::getMol() {
    return mol;
}