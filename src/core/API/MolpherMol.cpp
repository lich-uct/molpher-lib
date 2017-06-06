/*
 Copyright (c) 2016 Martin Šícho

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/BadFileException.h>
#include <boost/algorithm/string/predicate.hpp>

#include "data_structs/MolpherMol.hpp"
#include "MolpherMolImpl.hpp"
#include "core/misc/inout.h"

MolpherMol::MolpherMol(
    const std::string& string_repr
    , const std::string& formula
    , const std::string& parentSmile
    , const unsigned& oper
    , const double& dist
    , const double& distToClosestDecoy
    , const double& weight
    , const double& sascore)
    :
    pimpl(new MolpherMol::MolpherMolImpl(string_repr))
{
    pimpl->data.formula = formula;
    pimpl->data.parentSmile = parentSmile;
    pimpl->data.parentOper = oper;
    pimpl->data.distToTarget = dist;
    pimpl->data.molecularWeight = weight;
    pimpl->data.sascore = sascore;
}

MolpherMol::MolpherMol(const std::string& string_repr) : pimpl(new MolpherMol::MolpherMolImpl(string_repr)) {
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

std::shared_ptr<MolpherMol> MolpherMol::copy() const {
    return std::make_shared<MolpherMol>(*this);
}

MolpherMol::~MolpherMol() = default;

MolpherMol& MolpherMol::operator=(const MolpherMol& other) {
    pimpl = std::move(other.pimpl->copy());
}

// pimpl

MolpherMol::MolpherMolImpl::MolpherMolImpl() {
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const std::string& string_repr) {
    this->initialize(string_repr);
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMolData& data) : data(data) {
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMol::MolpherMolImpl& other) : data(other.data), tree(other.tree) {
    // no action
}

void MolpherMol::MolpherMolImpl::initialize(std::unique_ptr<RDKit::RWMol> mol) {
    try {
        RDKit::MolOps::Kekulize(*mol);
    } catch (const ValueErrorException &exc) {
        SynchCerr("Cannot kekulize input molecule.");
        throw exc;
    }

    data.SMILES = RDKit::MolToSmiles(*mol);
    data.formula = RDKit::Descriptors::calcMolFormula(*mol);
    rd_mol = std::move(mol);
}

void MolpherMol::MolpherMolImpl::initialize(const std::string &string_repr) {
    bool is_owned = (bool) tree;
    if (!is_owned) {
        if (string_repr.empty()) {
            SynchCerr("Creating a molecule with an empty SMILES string.");
            data.SMILES = "";
            return;
        }

        std::unique_ptr<RDKit::RWMol> mol;
        try {
            if (boost::algorithm::ends_with(string_repr, ".sdf")) {
                RDKit::ROMol* mol_ro = RDKit::SDMolSupplier(string_repr).next();
                mol.reset(new RDKit::RWMol(*mol_ro));
                delete mol_ro;
            } else if (string_repr.find("\n") != std::string::npos) {
                std::istream* istr = new std::istringstream(string_repr);
                RDKit::ROMol* mol_ro = RDKit::SDMolSupplier(istr).next();
                mol.reset(new RDKit::RWMol(*mol_ro));
                delete mol_ro;
            } else {
                mol.reset(RDKit::SmilesToMol(string_repr));
            }
        } catch (RDKit::SmilesParseException &exp) {
            SynchCerr("Error parsing supplied SMILES: \"" + string_repr + "\"");
            SynchCerr(exp.what());
            throw exp;
        }

        initialize(std::move(mol));
    } else {
        throw std::runtime_error("Molecule is already associated with a tree. "
            "The SMILES string cannot be changed at the moment.");
    }

//    SynchCout("Parsed molecule " + string_repr + " >> " + data.SMILES);
}

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
    pimpl->initialize(smiles);
}

void MolpherMol::setParentSMILES(const std::string& smiles) {
    bool is_owned = (bool) pimpl->tree;
    if (!is_owned) {
        pimpl->data.parentSmile = smiles;
    } else {
        throw std::runtime_error("Molecule is associated with a tree. "
            "SMILES of the parent cannot be modified.");
    }
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

void MolpherMol::removeFromTree() {
    auto tree = pimpl->tree;
    if (tree) {
        pimpl->tree.reset();
        auto smiles = getSMILES();
        if (tree->hasMol(smiles)) {
            tree->deleteSubtree(smiles, false);
        }
    }
}

