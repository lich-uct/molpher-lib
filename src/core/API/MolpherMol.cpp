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
#include <RDGeneral/BadFileException.h>
#include <boost/algorithm/string/predicate.hpp>
#include <core/chem/morphing/ReturnResults.hpp>
#include <core/misc/utils.hpp>

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
    , const double& sascore
    , const std::set<int>& fixed_atoms)
    :
    pimpl(new MolpherMol::MolpherMolImpl(string_repr))
{
    pimpl->data.formula = formula;
    pimpl->data.parentSmile = parentSmile;
    pimpl->data.parentOper = oper;
    pimpl->data.distToTarget = dist;
    pimpl->data.molecularWeight = weight;
    pimpl->data.sascore = sascore;
//    pimpl->fixed_atoms = fixed_atoms;
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

MolpherMol::MolpherMol(
        RDKit::ROMol* rd_mol
        , const std::string& formula
        , const std::string& parentSmile
        , const unsigned& oper
        , const double& dist
        , const double& distToClosestDecoy
        , const double& weight
        , const double& sascore
        , const std::set<int>& fixed_atoms
)
        : pimpl(new MolpherMol::MolpherMolImpl(
            *rd_mol
            , formula
            , parentSmile
            , oper
            , dist
            , distToClosestDecoy
            , weight
            , sascore
            , fixed_atoms
        ))
{
    // no action
}

MolpherMol::MolpherMol(RDKit::ROMol *rd_mol)
: pimpl(new MolpherMol::MolpherMolImpl(std::move(std::unique_ptr<RDKit::RWMol>(new RDKit::RWMol(*rd_mol)))))
{
    // no action
}

MolpherMol::MolpherMol(RDKit::RWMol *&rd_mol)
: pimpl(new MolpherMol::MolpherMolImpl(std::move(std::unique_ptr<RDKit::RWMol>(rd_mol))))
{
    rd_mol = nullptr;
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
    return *this;
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

MolpherMol::MolpherMolImpl::MolpherMolImpl(const MolpherMol::MolpherMolImpl& other)
:
data(other.data)
, tree(other.tree)
, rd_mol(new RDKit::RWMol(*other.rd_mol))
, atoms(other.atoms)
{
    // no action
}

MolpherMol::MolpherMolImpl::MolpherMolImpl(
        const RDKit::ROMol &rd_mol
        , const std::string& formula
        , const std::string& parentSmile
        , const unsigned& oper
        , const double& dist
        , const double& distToClosestDecoy
        , const double& weight
        , const double& sascore
        , const std::set<int>& fixed_atoms
) {
    data.formula = formula;
    data.parentSmile = parentSmile;
    data.parentOper = oper;
    data.distToTarget = dist;
    data.molecularWeight = weight;
    data.sascore = sascore;
//    this->fixed_atoms = fixed_atoms;

    std::unique_ptr<RDKit::RWMol> new_mol(new RDKit::RWMol(rd_mol));
    this->initialize(std::move(new_mol));
}

void MolpherMol::MolpherMolImpl::initialize(std::unique_ptr<RDKit::RWMol> mol) {
    try {
        if( !mol->getRingInfo()->isInitialized() ) {
            RDKit::MolOps::findSSSR(*mol);
        }
        RDKit::MolOps::sanitizeMol(*mol);
    } catch (const ValueErrorException &exc) {
        SynchCerr("Cannot kekulize input molecule.");
        throw exc;
    }

    RDKit::ROMol::AtomIterator iter;
    for (iter = mol->beginAtoms(); iter != mol->endAtoms(); iter++) {
        RDKit::Atom* atom = *iter;
        atoms.push_back(std::make_shared<MolpherAtom>(atom));
    }

    if (!mol->getPropList().empty()) {
        for (auto& item : parse_atom_locks(*mol)) {
            atoms[item.first]->setLockingMask(item.second);
        }
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

MolpherMol::MolpherMolImpl::MolpherMolImpl(std::unique_ptr<RDKit::RWMol> mol) {
    initialize(std::move(mol));
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

RDKit::RWMol* MolpherMol::asRDMol() const {
    return new RDKit::RWMol(*(pimpl->rd_mol));
}

void MolpherMol::lockAtom(int idx, int mask) {
    getAtom(idx)->setLockingMask(mask);
}

std::shared_ptr<MolpherAtom> MolpherMol::getAtom(int idx) const {
    return pimpl->atoms[idx];
}

const std::vector<std::shared_ptr<MolpherAtom>> &MolpherMol::getAtoms() const {
    return pimpl->atoms;
}

int MolpherMol::getAtomCount() const {
	return pimpl->atoms.size();
}
