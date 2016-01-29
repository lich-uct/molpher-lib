
#include "data_structs/MolpherMol.hpp"

MolpherMol::MolpherMol(std::shared_ptr<MolpherMolImpl> pimpl) : pimpl(pimpl) {
    // no action
}



//MolpherMol::MolpherMol() : mol(new MolpherMolecule()), selfAllocated(true) {
//    // no action
//}
//
//
//MolpherMol::MolpherMol(const MolpherMol& other) : mol(nullptr), selfAllocated(true) {
//    if (other.selfAllocated) {
//        MolpherMolecule* new_mol = new MolpherMolecule();
//        *new_mol = *(other.mol);
//        this->mol = new_mol;
//    } else {
//        selfAllocated = false;
//        this->mol = other.mol;
//    }
//}
//
//MolpherMol::MolpherMol(MolpherMolecule& mol) : mol(&mol), selfAllocated(false) {
//    // no action
//}
//
//MolpherMol::MolpherMol(MolpherMolecule& mol, bool copy) : mol(nullptr), selfAllocated(false) {
//    if (copy) {
//        MolpherMolecule* new_mol = new MolpherMolecule();
//        *new_mol = mol;
//        this->mol = new_mol;
//        selfAllocated = true;
//    } else {
//        this->mol = &mol;
//        selfAllocated = false;
//    }
//}
//
//MolpherMol::~MolpherMol() {
//    if (selfAllocated) delete mol;
//}
//
//MolpherMol& MolpherMol::operator=(const MolpherMol& other) {
//    if (selfAllocated) {
//        MolpherMolecule* new_mol = new MolpherMolecule();
//        *new_mol = *(other.mol);
//        delete this->mol;
//        this->mol = nullptr;
//        this->mol = new_mol;
//        this->selfAllocated = true;
//    } else {
//        this->mol = other.mol;
//    }
//    return *this;
//}
//
//MolpherMolecule& MolpherMol::fetchMolpherMolecule() const {
//    return *mol;
//}
//
//bool MolpherMol::isBound() const {
//    return !selfAllocated;
//}
//
//MolpherMol* MolpherMol::copy() const {
//    MolpherMolecule* new_mol = new MolpherMolecule();
//    *new_mol = *mol;
//    MolpherMol* ret = new MolpherMol(*new_mol);
//    ret->selfAllocated = true;
//    return ret;
//}
//
//std::string MolpherMol::getSMILES() {
//    return mol->smile;
//}
//
//double MolpherMol::getDistToTarget() {
//    return mol->distToTarget;
//}
//
//std::string MolpherMol::getParentSMILES() {
//    return mol->parentSmile;
//}
//
//const std::set<std::string>& MolpherMol::getDescendants() {
//    return mol->descendants;
//}
//
//const std::set<std::string>& MolpherMol::getHistoricDescendants() {
//    return mol->historicDescendants;
//}
//
//unsigned int MolpherMol::getItersWithoutDistImprovement() {
//    return (unsigned int) mol->itersWithoutDistImprovement;
//}
//
//double MolpherMol::getMolecularWeight() {
//    return mol->molecularWeight;
//}
//
//double MolpherMol::getSAScore() {
//    return mol->sascore;
//}
//
//void MolpherMol::setDistToTarget(double dist) {
//    mol->distToTarget = dist;
//}
//
//void MolpherMol::setSAScore(double dist) {
//    mol->sascore = dist;
//}

