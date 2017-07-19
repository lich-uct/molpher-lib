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
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>
#include <core/chem/morphing/ReturnResults.hpp>

#include "MolpherMolImpl.hpp"
#include "core/misc/inout.h"
#include "core/API/operations/GenerateMorphsOperImpl.hpp"
#include "core/chem/morphing/CalculateMorphs.hpp"
#include "core/chem/morphing/CalculateDistances.hpp"

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
    pimpl->fixed_atoms = fixed_atoms;
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
        RDKit::RWMol* rd_mol
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

MolpherMol::MolpherMolImpl::MolpherMolImpl(
        const RDKit::RWMol &rd_mol
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
    this->fixed_atoms = fixed_atoms;

    std::unique_ptr<RDKit::RWMol> new_mol(new RDKit::RWMol(rd_mol));
    this->initialize(std::move(new_mol));
}

void MolpherMol::MolpherMolImpl::initialize(std::unique_ptr<RDKit::RWMol> mol) {
    try {
        if( !mol->getRingInfo()->isInitialized() ) {
            RDKit::MolOps::findSSSR(*mol);
        }
        RDKit::MolOps::Kekulize(*mol);
    } catch (const ValueErrorException &exc) {
        SynchCerr("Cannot kekulize input molecule.");
        throw exc;
    }

    RDKit::STR_VECT prop_names = mol->getPropList();
    if (std::find(prop_names.begin(), prop_names.end(), "MOLPHER_FIXED") != prop_names.end()) {
        std::string fixed_positions = mol->getProp<std::string>("MOLPHER_FIXED");
        std::vector<std::string> indices;
        split(fixed_positions, ',', std::back_inserter(indices));
        for (auto idx : indices) {
            fixed_atoms.insert(std::stoi(idx) - 1);
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

std::unique_ptr<ConcurrentMolVector>
MolpherMol::MolpherMolImpl::morph(const std::vector<ChemOperSelector> &operators, int cntMorphs, int threadCnt,
                                  FingerprintSelector fingerprintSelector, SimCoeffSelector simCoeffSelector,
                                  const MolpherMol &target) {
    tbb::task_group_context tbbCtx;
    tbb::task_scheduler_init scheduler;
    scheduler.terminate();
    scheduler.initialize(threadCnt);

    std::unique_ptr<ConcurrentMolVector> candidates(new ConcurrentMolVector());
    candidates->reserve(candidates->size() + cntMorphs);
    CollectMorphs collectMorphs(*candidates, false);

    RDKit::ROMol* targetMol = target.pimpl->rd_mol.get();
    SimCoefCalculator scCalc(simCoeffSelector , fingerprintSelector, rd_mol.get(), targetMol);
    Fingerprint* targetFp = scCalc.GetFingerprint(targetMol);

    std::vector<MorphingStrategy *> strategies;
    InitStrategies(operators, strategies);

    RDKit::RWMol **newMols = new RDKit::RWMol *[cntMorphs];
    std::memset(newMols, 0, sizeof(RDKit::RWMol *) * cntMorphs);
    ChemOperSelector *opers = new ChemOperSelector [cntMorphs];
    std::string *smiles = new std::string [cntMorphs];
    std::string *formulas = new std::string [cntMorphs];
    double *weights = new double [cntMorphs];
    double *sascores = new double [cntMorphs]; // added for SAScore
    double *distToTarget = new double [cntMorphs];

    // compute new morphs and smiles
    if (!tbbCtx.is_group_execution_cancelled()) {
        tbb::atomic<unsigned int> kekulizeFailureCount;
        tbb::atomic<unsigned int> sanitizeFailureCount;
        tbb::atomic<unsigned int> morphingFailureCount;
        kekulizeFailureCount = 0;
        sanitizeFailureCount = 0;
        morphingFailureCount = 0;
        try {
            MorphingData data(*rd_mol, *targetMol, operators);
            CalculateMorphs calculateMorphs(
                    data, strategies, opers, newMols, smiles, formulas, weights, sascores,
                    kekulizeFailureCount, sanitizeFailureCount, morphingFailureCount);

            tbb::parallel_for(tbb::blocked_range<int>(0, cntMorphs),
                              calculateMorphs, tbb::auto_partitioner(), tbbCtx);
        } catch (const std::exception &exc) {
//            REPORT_RECOVERY("Recovered from morphing data construction failure.");
        }
        if (kekulizeFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << kekulizeFailureCount << " kekulization failures.";
//            REPORT_RECOVERY(report.str());
        }
        if (sanitizeFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << sanitizeFailureCount << " sanitization failures.";
//            REPORT_RECOVERY(report.str());
        }
        if (morphingFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << morphingFailureCount << " morphing failures.";
//            REPORT_RECOVERY(report.str());
        }
    }

    // compute distances
    // we need to announce the decoy which we want to use
    if (!tbbCtx.is_group_execution_cancelled()) {
        CalculateDistances calculateDistances(newMols, scCalc, targetFp, distToTarget);
        tbb::parallel_for(tbb::blocked_range<int>(0, cntMorphs),
                          calculateDistances, tbb::auto_partitioner(), tbbCtx);
    }

    // return results
    if (!tbbCtx.is_group_execution_cancelled()) {
        std::string parent = data.SMILES;
        ReturnResults returnResults(
                newMols, smiles, formulas, parent, opers, weights, sascores,
                distToTarget, &collectMorphs, CollectMorphs::MorphCollector);
        tbb::parallel_for(tbb::blocked_range<int>(0, cntMorphs),
                          returnResults, tbb::auto_partitioner(), tbbCtx);
    }

    // clean up
    for (int i = 0; i < cntMorphs; ++i) {
        delete newMols[i];
    }
    delete[] newMols;
    delete[] opers;
    delete[] smiles;
    delete[] formulas;
    delete[] weights;
    delete[] sascores;

    delete targetFp;

    for (int i = 0; i < strategies.size(); ++i) {
        delete strategies[i];
    }

    return std::move(candidates);
}

std::unique_ptr<ConcurrentMolVector>
MolpherMol::MolpherMolImpl::morph(const std::vector<ChemOperSelector> &operators, int cntMorphs, int threadCnt) {
    tbb::task_group_context tbbCtx;
    tbb::task_scheduler_init scheduler;
    scheduler.terminate();
    scheduler.initialize(threadCnt);

    std::unique_ptr<ConcurrentMolVector> candidates(new ConcurrentMolVector());
    candidates->reserve(candidates->size() + cntMorphs);
    CollectMorphs collectMorphs(*candidates, false);

    std::vector<MorphingStrategy *> strategies;
    InitStrategies(operators, strategies);

    RDKit::RWMol **newMols = new RDKit::RWMol *[cntMorphs];
    ChemOperSelector *opers = new ChemOperSelector [cntMorphs];
    std::memset(newMols, 0, sizeof(RDKit::RWMol *) * cntMorphs);
    std::string *smiles = new std::string [cntMorphs];
    std::string *formulas = new std::string [cntMorphs];
    double *weights = new double [cntMorphs];
    double *sascores = new double [cntMorphs]; // added for SAScore

    // compute new morphs and smiles
    MorphingData morphing_data(*rd_mol, operators, fixed_atoms);
    if (!tbbCtx.is_group_execution_cancelled()) {
        tbb::atomic<unsigned int> kekulizeFailureCount;
        tbb::atomic<unsigned int> sanitizeFailureCount;
        tbb::atomic<unsigned int> morphingFailureCount;
        kekulizeFailureCount = 0;
        sanitizeFailureCount = 0;
        morphingFailureCount = 0;
        try {
            CalculateMorphs calculateMorphs(
                    morphing_data, strategies, opers, newMols, smiles, formulas, weights, sascores,
                    kekulizeFailureCount, sanitizeFailureCount, morphingFailureCount);

            tbb::parallel_for(tbb::blocked_range<int>(0, cntMorphs),
                              calculateMorphs, tbb::auto_partitioner(), tbbCtx);
        } catch (const std::exception &exc) {
//            REPORT_RECOVERY("Recovered from morphing data construction failure.");
        }
        if (kekulizeFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << kekulizeFailureCount << " kekulization failures.";
//            REPORT_RECOVERY(report.str());
        }
        if (sanitizeFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << sanitizeFailureCount << " sanitization failures.";
//            REPORT_RECOVERY(report.str());
        }
        if (morphingFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << morphingFailureCount << " morphing failures.";
//            REPORT_RECOVERY(report.str());
        }
    }

    // return results
    if (!tbbCtx.is_group_execution_cancelled()) {
        std::string parent = data.SMILES;
        ReturnResults returnResults(
                newMols, smiles, formulas, parent, opers, weights, sascores,
                fixed_atoms,
                morphing_data.removed_atoms,
                &collectMorphs, CollectMorphs::MorphCollector);
        tbb::parallel_for(tbb::blocked_range<int>(0, cntMorphs),
                          returnResults, tbb::auto_partitioner(), tbbCtx);
    }

    // clean up
    for (int i = 0; i < cntMorphs; ++i) {
        delete newMols[i];
    }
    delete[] newMols;
    delete[] opers;
    delete[] smiles;
    delete[] formulas;
    delete[] weights;
    delete[] sascores;

    for (int i = 0; i < strategies.size(); ++i) {
        delete strategies[i];
    }

    return std::move(candidates);
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

std::vector<std::shared_ptr<MolpherMol> >
MolpherMol::morph(const std::vector<ChemOperSelector> &operators, int cntMorphs, int threadCnt,
				  FingerprintSelector fingerprintSelector, SimCoeffSelector simCoeffSelector,
				  const MolpherMol &target) {
    MolVector ret;
    std::unique_ptr<ConcurrentMolVector> candidates(pimpl->morph(operators, cntMorphs, threadCnt, fingerprintSelector, simCoeffSelector, target));
    for (auto morph : *candidates) {
        ret.push_back(morph);
    }
    return ret;
}

std::vector<std::shared_ptr<MolpherMol> >
MolpherMol::morph(const std::vector<ChemOperSelector> &operators, int cntMorphs, int threadCnt) {
    MolVector ret;
    std::unique_ptr<ConcurrentMolVector> candidates(pimpl->morph(operators, cntMorphs, threadCnt));
    for (auto morph : *candidates) {
        ret.push_back(morph);
    }
    return ret;
}
