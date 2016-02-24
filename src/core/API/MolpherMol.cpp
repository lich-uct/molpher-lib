
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/BadFileException.h>

#include "data_structs/MolpherMol.hpp"
#include "MolpherMolImpl.hpp"
#include "data_structs/ExplorationTree.hpp"
#include "core/misc/inout.h"

//MolpherMol::MolpherMol(std::shared_ptr<MolpherMolImpl> pimpl) : pimpl(pimpl) {
//    // no action
//}

MolpherMol::MolpherMol(
    const std::string& smiles
    , const std::string& formula
    , const std::string& parentSmile
    , const unsigned& oper
    , const double& dist
    , const double& distToClosestDecoy
    , const double& weight
    , const double& sascore)
    :
    pimpl(new MolpherMol::MolpherMolImpl(smiles))
{
    pimpl->data.formula = formula;
    pimpl->data.parentSmile = parentSmile;
    pimpl->data.parentOper = oper;
    pimpl->data.distToTarget = dist;
    pimpl->data.molecularWeight = weight;
    pimpl->data.sascore = sascore;
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

std::shared_ptr<ExplorationTree> MolpherMol::getTree() {
    return pimpl->tree;
}

std::unique_ptr<MolpherMol> MolpherMol::copy() const {
    return std::unique_ptr<MolpherMol>(new MolpherMol(*this));
}

MolpherMol::~MolpherMol() = default;

MolpherMol& MolpherMol::operator=(const MolpherMol& other) {
    pimpl = std::move(other.pimpl->copy());
}

// pimpl

MolpherMol::MolpherMolImpl::MolpherMolImpl() {
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const std::string& smiles) {
    this->setSMILES(smiles);
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMolData& data) : data(data) {
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMol::MolpherMolImpl& other) : data(other.data), tree(other.tree) {
    // no action
}

void MolpherMol::MolpherMolImpl::setSMILES(const std::string& smiles) {
    RDKit::RWMol* mol = nullptr;
    try {
        if (smiles.empty()) {
            SynchCerr("Creating a molecule with an empty SMILES string.");
            data.SMILES = "";
            return;
        }
        mol = RDKit::SmilesToMol(smiles);
    } catch (RDKit::SmilesParseException &exp) {
        SynchCerr("Error parsing supplied SMILES: \"" + smiles + "\"");
        SynchCerr(exp.what());
        throw exp;
    }
    
    try {
        RDKit::MolOps::Kekulize(*mol);
    } catch (const ValueErrorException &exc) {
        SynchCerr("Cannot kekulize input molecule.");
        throw exc;
    }

    data.SMILES = RDKit::MolToSmiles(*mol);
    data.formula = RDKit::Descriptors::calcMolFormula(*mol);

    SynchCout("Parsed molecule " + smiles + " >> " + data.SMILES);
}

//MolpherMolData MolpherMol::MolpherMolImpl::asData() const {
//    return data;
//}

const std::string& MolpherMol::getSMILES() const {
    return pimpl->data.SMILES;
}

void MolpherMol::setDistToTarget(double dist) {
    pimpl->data.distToTarget = dist;
}

double MolpherMol::getDistToTarget() const {
    return pimpl->data.distToTarget;
}

std::unique_ptr<MolpherMol::MolpherMolImpl> MolpherMol::MolpherMolImpl::copy() const {
    return std::unique_ptr<MolpherMol::MolpherMolImpl>(new MolpherMol::MolpherMolImpl(*this));
}

void MolpherMol::addToDescendants(const std::string& smiles) {
    pimpl->data.descendants.insert(smiles);
}

void MolpherMol::addToHistoricDescendants(const std::string& smiles) {
    pimpl->data.historicDescendants.insert(smiles);
}

void MolpherMol::decreaseItersWithoutDistImprovement() {
    pimpl->data.gensWithoutDistImprovement--;
}

const std::set<std::string>& MolpherMol::getDescendants() const {
    return pimpl->data.descendants;
}

const std::string& MolpherMol::getFormula() const {
    return pimpl->data.formula;
}

const std::set<std::string>& MolpherMol::getHistoricDescendants() const {
    return pimpl->data.historicDescendants;
}

unsigned int MolpherMol::getItersWithoutDistImprovement() const {
    return pimpl->data.gensWithoutDistImprovement;
}

double MolpherMol::getMolecularWeight() const {
    return pimpl->data.molecularWeight;
}

int MolpherMol::getParentOper() const {
    return pimpl->data.parentOper;
}

const std::string& MolpherMol::getParentSMILES() const {
    return pimpl->data.parentSmile;
}

double MolpherMol::getSAScore() const {
    return pimpl->data.sascore;
}

void MolpherMol::increaseItersWithoutDistImprovement() {
    pimpl->data.gensWithoutDistImprovement++;
}

bool MolpherMol::isBoundToTree() const {
    return (bool) pimpl->tree;
}

bool MolpherMol::isValid() const {
    return pimpl->data.isValid();
}

void MolpherMol::removeFromDescendants(const std::string& smiles) {
    auto& descs = pimpl->data.descendants;
    if (descs.find(smiles) != descs.end()) {
        descs.erase(smiles);
    } else {
        SynchCerr("Molecule (" + smiles + ") not found in descendants. No changes made...");
    }
}

void MolpherMol::removeFromHistoricDescendants(const std::string& smiles) {
    auto& descs = pimpl->data.historicDescendants;
    if (descs.find(smiles) != descs.end()) {
        descs.erase(smiles);
    } else {
        SynchCerr("Molecule (" + smiles + ") not found in historic descendants. No changes made...");
    }
}

void MolpherMol::setDescendants(const std::set<std::string>& new_set) {
    pimpl->data.descendants = new_set;
}

void MolpherMol::setHistoricDescendants(const std::set<std::string>& new_set) {
    pimpl->data.historicDescendants = new_set;
}

void MolpherMol::setItersWithoutDistImprovement(unsigned int count) {
    pimpl->data.gensWithoutDistImprovement = count;
}

void MolpherMol::setSAScore(double score) {
    pimpl->data.sascore = score;
}

void MolpherMol::setSMILES(const std::string& smiles) {
    pimpl->setSMILES(smiles);
}

void MolpherMol::setParentSMILES(const std::string& smiles) {
    pimpl->data.parentSmile = smiles;
}

void MolpherMol::setOwner(std::shared_ptr<ExplorationTree> tree) {
    bool is_owned = (bool) pimpl->tree;
    if (!is_owned) {
        pimpl->tree = tree;
    } else {
        throw std::runtime_error("Molecule is already associated with a tree. "
            "Molecules cannot be assigned more than once.");
    }
}

