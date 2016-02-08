
#include "data_structs/MolpherMol.hpp"
#include "MolpherMolImpl.hpp"
#include "data_structs/ExplorationTree.hpp"

//MolpherMol::MolpherMol(std::shared_ptr<MolpherMolImpl> pimpl) : pimpl(pimpl) {
//    // no action
//}

MolpherMol::MolpherMol(std::string& smiles, std::string& formula, std::string& parentSmile, ChemOperSelector* opers, double dist, double distToClosestDecoy, double weight, double sascore) {
    // TODO: implement
}

MolpherMol::MolpherMol(const std::string& smiles) : pimpl(new MolpherMol::MolpherMolImpl(smiles)) {
    // no action
}

MolpherMol::MolpherMol() : pimpl(new MolpherMol::MolpherMolImpl()) {
    // no action
}

MolpherMol::MolpherMol(const MolpherMol& other) : pimpl(std::move(other.pimpl->copy())) {
    // no action
}

std::string MolpherMol::getSMILES() {
    return pimpl->getSMILES();
}

void MolpherMol::setDistToTarget(double dist) {
    pimpl->setDistToTarget(dist);
}

double MolpherMol::getDistToTarget() {
    return pimpl->getDistToTarget();
}

std::shared_ptr<ExplorationTree> MolpherMol::getTree() {
    return pimpl->getTree();
}

std::unique_ptr<MolpherMol> MolpherMol::copy() const {
    return std::unique_ptr<MolpherMol>(new MolpherMol(*this));
}

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

MolpherMol::~MolpherMol() = default;

MolpherMol& MolpherMol::operator=(const MolpherMol& other) {
    pimpl = std::move(other.pimpl->copy());
}

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

MolpherMol::MolpherMolImpl::MolpherMolImpl() {
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const std::string& smiles) {
    data.SMILES = smiles;
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMolData& data) : data(data) {
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMol::MolpherMolImpl& other) : data(other.data), tree(other.tree) {
    // no action
}

std::string MolpherMol::MolpherMolImpl::getSMILES() const {
    return data.SMILES;
}

void MolpherMol::MolpherMolImpl::setDistToTarget(double dist) {
    data.distToTarget = dist;
}

double MolpherMol::MolpherMolImpl::getDistToTarget() const {
    return data.distToTarget;
}

std::shared_ptr<ExplorationTree> MolpherMol::MolpherMolImpl::getTree() {
    return tree;
}


MolpherMolData MolpherMol::MolpherMolImpl::asData() const {
    return data;
}

std::unique_ptr<MolpherMol::MolpherMolImpl> MolpherMol::MolpherMolImpl::copy() const {
    return std::unique_ptr<MolpherMol::MolpherMolImpl>(new MolpherMol::MolpherMolImpl(*this));
}
