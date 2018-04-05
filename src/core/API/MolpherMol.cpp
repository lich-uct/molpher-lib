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
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/BadFileException.h>
#include <boost/algorithm/string/predicate.hpp>
#include <core/misc/utils.hpp>
#include <core/misc/SAScore.h>

#include "MolpherMolImpl.hpp"
#include "core/misc/inout.h"

MolpherMol::MolpherMol(
        const std::string& string_repr
        , const std::string& parentSmile
        , const std::string& parentOper
        , const double& dist
        , const double& sascore
)
    :
    pimpl(new MolpherMol::MolpherMolImpl(string_repr))
{
    pimpl->data.parentSmile = parentSmile;
    pimpl->data.parentOper = parentOper;
    pimpl->data.distToTarget = dist;
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

MolpherMol::MolpherMol(RDKit::ROMol *rd_mol)
: pimpl(new MolpherMol::MolpherMolImpl(std::move(std::unique_ptr<RDKit::RWMol>(new RDKit::RWMol(*rd_mol)))))
{
    // no action
}

MolpherMol::MolpherMol(RDKit::RWMol* &rd_mol)
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

void MolpherMol::MolpherMolImpl::initialize(std::unique_ptr<RDKit::RWMol> mol) {
    // sanitize
    unsigned int failed = 0;
    try {
        RDKit::MolOps::sanitizeMol(
                *mol
                , failed
                , RDKit::MolOps::SANITIZE_KEKULIZE
                | RDKit::MolOps::SANITIZE_PROPERTIES
                | RDKit::MolOps::SANITIZE_SYMMRINGS
                | RDKit::MolOps::SANITIZE_CLEANUP
                | RDKit::MolOps::SANITIZE_FINDRADICALS
                | RDKit::MolOps::SANITIZE_SETHYBRIDIZATION
                | RDKit::MolOps::SANITIZE_CLEANUPCHIRALITY
                | RDKit::MolOps::SANITIZE_ADJUSTHS
				| RDKit::MolOps::SANITIZE_SETCONJUGATION
        );
    } catch (const RDKit::MolSanitizeException &exc) {
        SynchCerr("Molecule failed to initialize due to sanitization errors.", "ERROR:");
        throw exc;
    } catch (const ValueErrorException &exc) {
        SynchCerr("Cannot kekulize input molecule.");
        throw exc;
    } catch (const std::exception &exc) {
        SynchCerr("Molecule failed to initialize due to an exception: " + std::string(exc.what()), "ERROR:");
        throw exc;
    }

    rd_mol = std::move(mol); // take ownership of the instance

    // transfer locks and calculate valences for atoms
	RDKit::ROMol::AtomIterator iter;
    for (iter = rd_mol->beginAtoms(); iter != rd_mol->endAtoms(); iter++) {
        RDKit::Atom* atom = *iter;
		atom->calcExplicitValence();
		atom->calcImplicitValence();
        atoms.push_back(std::make_shared<MolpherAtom>(atom));
    }
    if (!rd_mol->getPropList().empty()) {
        for (auto& item : parse_atom_locks(*rd_mol)) {
            lockAtom(item.first, item.second);
        }
    }

    data.SMILES = RDKit::MolToSmiles(*rd_mol, false, true, -1, true, false, false);
    data.formula = RDKit::Descriptors::calcMolFormula(*rd_mol);
    data.molecularWeight = RDKit::Descriptors::calcAMW(*rd_mol, true);
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

void MolpherMol::MolpherMolImpl::lockAtom(int idx, int mask) {
    std::shared_ptr<MolpherAtom> atom = getAtom(idx);
    atom->setLockingMask(mask);
    mask = atom->getLockingMask();

    if (mask & (MolpherAtom::KEEP_NEIGHBORS | MolpherAtom::KEEP_BONDS)) {
        for (auto ngh : getNeighbors(idx)) {
            if (mask & MolpherAtom::KEEP_NEIGHBORS) {
                ngh->setLockingMask(ngh->getLockingMask() | MolpherAtom::NO_MUTATION | MolpherAtom::NO_REMOVAL);
            }
            if (mask & MolpherAtom::KEEP_BONDS) {
                ngh->setLockingMask(ngh->getLockingMask() | MolpherAtom::NO_REMOVAL);
            }
        }
    }
}

std::shared_ptr<MolpherAtom> MolpherMol::MolpherMolImpl::getAtom(int idx) {
    if (idx < 0 || idx >= atoms.size()) {
        std::runtime_error("No such atom. Index out of range: " + std::to_string(idx));
    }
    return atoms[idx];
}

std::vector<std::shared_ptr<MolpherAtom>> MolpherMol::MolpherMolImpl::getNeighbors(int idx) {
    if (idx < 0 || idx >= atoms.size()) {
        std::runtime_error("No such atom. Index out of range: " + std::to_string(idx));
    }

    RDKit::ROMol& mol = *rd_mol;
    std::vector<std::shared_ptr<MolpherAtom>> neighbors;
    RDKit::ROMol::ADJ_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomNeighbors(mol.getAtomWithIdx(idx));
    while (beg != end) {
        RDKit::Atom *neighbor = mol[*beg].get();
        neighbors.push_back(getAtom(neighbor->getIdx()));
        ++beg;
    }
    return neighbors;
}

std::string MolpherMol::MolpherMolImpl::asMolBlock() const {
    return RDKit::MolToMolBlock(*rd_mol, true, -1, false, false);
}

std::shared_ptr<MolpherMol> MolpherMol::MolpherMolImpl::fromMolBlock(const std::string &mol_block) {
    RDKit::RWMol* mol = RDKit::MolBlockToMol(mol_block, false, false, false);
    return std::shared_ptr<MolpherMol>(new MolpherMol(mol));
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

const std::string& MolpherMol::getParentOper() const {
    return pimpl->data.parentOper;
}

const std::string& MolpherMol::getParentSMILES() const {
    return pimpl->data.parentSmile;
}

double MolpherMol::getSAScore() const {
    if (pimpl->data.sascore == 0) {
        // calculate SAScore
        // FIXME: sometimes we get scores that are too high if we use the rd_mol instance directly
        auto rd_mol_dummy = RDKit::SmilesToMol(getSMILES());
        pimpl->data.sascore = SAScore::getInstance()->getScore(*rd_mol_dummy);
        delete rd_mol_dummy;
    }
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
    return new RDKit::RWMol(*(pimpl->rd_mol), true, -1);
}

void MolpherMol::lockAtom(int idx, int mask) {
    pimpl->lockAtom(idx, mask);
}

std::shared_ptr<MolpherAtom> MolpherMol::getAtom(int idx) const {
    return pimpl->getAtom(idx);
}

const std::vector<std::shared_ptr<MolpherAtom>> &MolpherMol::getAtoms() const {
    return pimpl->atoms;
}

int MolpherMol::getAtomCount() const {
	return pimpl->atoms.size();
}

const std::vector<std::shared_ptr<MolpherAtom>> MolpherMol::getNeighbors(int idx) const {
    return pimpl->getNeighbors(idx);
}

std::string MolpherMol::asMolBlock() const {
	return pimpl->asMolBlock();
}

std::shared_ptr<MolpherMol> MolpherMol::fromMolBlock(const std::string &mol_block) {
    return MolpherMolImpl::fromMolBlock(mol_block);
}

void MolpherMol::setParentOper(const std::string& description) {
    pimpl->data.parentOper = description;
}