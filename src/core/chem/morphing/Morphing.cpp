/*
 Copyright (c) 2012 Peter Szepe

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

#include <cstring>
#include <string>
#include <sstream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <tbb/atomic.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

//#include "main.hpp"
#include "core/misc/inout.h"
#include "core/chem/fingerprintStrategy/FingerprintStrategy.h"
#include "core/chem/simCoefStrategy/SimCoefStrategy.h"
#include "MorphingData.h"
#include "core/chem/morphingStrategy/MorphingStrategy.h"
#include "core/chem/morphing/MorphingFtors.hpp"
#include "core/chem/morphing/Morphing.hpp"

// TODO: merge into one header file ?
#include "core/chem/morphingStrategy/OpAddAtom.hpp"
#include "core/chem/morphingStrategy/OpAddBond.hpp"
#include "core/chem/morphingStrategy/OpBondContraction.hpp"
#include "core/chem/morphingStrategy/OpBondReroute.hpp"
#include "core/chem/morphingStrategy/OpInterlayAtom.hpp"
#include "core/chem/morphingStrategy/OpMutateAtom.hpp"
#include "core/chem/morphingStrategy/OpRemoveAtom.hpp"
#include "core/chem/morphingStrategy/OpRemoveBond.hpp"

#if MORPHING_REPORTING == 1
#define REPORT_RECOVERY(x) SynchCout((x))
#else
#define REPORT_RECOVERY(x)
#endif

static void InitStrategies(
    std::vector<ChemOperSelector> &chemOperSelectors,
    std::vector<MorphingStrategy *> &strategies)
{
    for (int i = 0; i < chemOperSelectors.size(); ++i) {
        switch (chemOperSelectors[i]) {
        case OP_ADD_ATOM:
            strategies.push_back(new OpAddAtom());
            break;
        case OP_REMOVE_ATOM:
            strategies.push_back(new OpRemoveAtom());
            break;
        case OP_ADD_BOND:
            strategies.push_back(new OpAddBond());
            break;
        case OP_REMOVE_BOND:
            strategies.push_back(new OpRemoveBond());
            break;
        case OP_MUTATE_ATOM:
            strategies.push_back(new OpMutateAtom());
            break;
        case OP_INTERLAY_ATOM:
            strategies.push_back(new OpInterlayAtom());
            break;
        case OP_BOND_REROUTE:
            strategies.push_back(new OpBondReroute());
            break;
        case OP_BOND_CONTRACTION:
            strategies.push_back(new OpBondContraction());
            break;
        default:
            break;
        }
    }
}

void GenerateMorphs(
    MolpherMol &candidate,
    unsigned int morphAttempts,
    FingerprintSelector fingerprintSelector,
    SimCoeffSelector simCoeffSelector,
    std::vector<ChemOperSelector> &chemOperSelectors,
    MolpherMol&target,
    std::vector<MolpherMol> &decoys,
    tbb::task_group_context &tbbCtx ,
    void *callerState,
    void (*deliver)(std::shared_ptr<MolpherMol>, void *)
) {
    RDKit::RWMol *mol = NULL;
    try {
        mol = RDKit::SmilesToMol(candidate.getSMILES());
        if (mol) {
            RDKit::MolOps::Kekulize(*mol);
        } else {
            throw ValueErrorException("");
        }
    } catch (const ValueErrorException &exc) {
        delete mol;
        return;
    }

    RDKit::RWMol *targetMol = NULL;
    try {
        targetMol = RDKit::SmilesToMol(target.getSMILES());
        if (targetMol) {
            RDKit::MolOps::Kekulize(*targetMol);
        } else {
            throw ValueErrorException("");
        }
    } catch (const ValueErrorException &exc) {
        delete targetMol;
        delete mol;
        return;
    }

    SimCoefCalculator scCalc(simCoeffSelector , fingerprintSelector, mol, targetMol);

    Fingerprint *targetFp = scCalc.GetFingerprint(targetMol);

    std::vector<Fingerprint *> decoysFp;
    decoysFp.reserve(decoys.size());
    RDKit::RWMol *decoyMol = NULL;
    try {
        for (int i = 0; i < decoys.size(); ++i) {
            RDKit::RWMol *decoyMol = RDKit::SmilesToMol(decoys[i].getSMILES());
            if (decoyMol) {
                RDKit::MolOps::Kekulize(*decoyMol);
                decoysFp.push_back(scCalc.GetFingerprint(decoyMol));
                delete decoyMol;
                decoyMol = NULL;
            } else {
                throw ValueErrorException("");
            }
        }
    } catch (const ValueErrorException &exc) {
        REPORT_RECOVERY("Recovered from decoy kekulization failure.");
        for (int i = 0; i < decoysFp.size(); ++i) {
            delete decoysFp[i];
        }
        delete decoyMol;
        delete targetMol;
        delete mol;
        return;
    }
    
    std::vector<MorphingStrategy *> strategies;
    InitStrategies(chemOperSelectors, strategies);

    RDKit::RWMol **newMols = new RDKit::RWMol *[morphAttempts];
    std::memset(newMols, 0, sizeof(RDKit::RWMol *) * morphAttempts);
    ChemOperSelector *opers = new ChemOperSelector [morphAttempts];
    std::string *smiles = new std::string [morphAttempts];
    std::string *formulas = new std::string [morphAttempts];
    double *weights = new double [morphAttempts];
    double *sascores = new double [morphAttempts]; // added for SAScore
    double *distToTarget = new double [morphAttempts];
    double *distToClosestDecoy = new double [morphAttempts];
                
    // compute new morphs and smiles
    if (!tbbCtx.is_group_execution_cancelled()) {
        tbb::atomic<unsigned int> kekulizeFailureCount;
        tbb::atomic<unsigned int> sanitizeFailureCount;
        tbb::atomic<unsigned int> morphingFailureCount;
        kekulizeFailureCount = 0;
        sanitizeFailureCount = 0;
        morphingFailureCount = 0;
        try {
			// TODO: incorporate this under the MolpherMol class -> MolpherMol instances generate their own morphs

			MorphingData data(*mol, *targetMol, chemOperSelectors);
            CalculateMorphs calculateMorphs(
                data, strategies, opers, newMols, smiles, formulas, weights, sascores,
                kekulizeFailureCount, sanitizeFailureCount, morphingFailureCount);
            
            tbb::parallel_for(tbb::blocked_range<int>(0, morphAttempts),
                calculateMorphs, tbb::auto_partitioner(), tbbCtx);
        } catch (const std::exception &exc) {
            REPORT_RECOVERY("Recovered from morphing data construction failure.");
        }
        if (kekulizeFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << kekulizeFailureCount << " kekulization failures.";
            REPORT_RECOVERY(report.str());
        }
        if (sanitizeFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << sanitizeFailureCount << " sanitization failures.";
            REPORT_RECOVERY(report.str());
        }
        if (morphingFailureCount > 0) {
            std::stringstream report;
            report << "Recovered from " << morphingFailureCount << " morphing failures.";
            REPORT_RECOVERY(report.str());
        }
    }
    
    // compute distances
    // we need to announce the decoy which we want to use    
    if (!tbbCtx.is_group_execution_cancelled()) {
        CalculateDistances calculateDistances(newMols, scCalc, targetFp,
            decoysFp, distToTarget, distToClosestDecoy, 0/*candidate.nextDecoy*/);
        tbb::parallel_for(tbb::blocked_range<int>(0, morphAttempts),
            calculateDistances, tbb::auto_partitioner(), tbbCtx);
    }

    // return results
    if (!tbbCtx.is_group_execution_cancelled()) {
        std::string parent = candidate.getSMILES();
        ReturnResults returnResults(
            newMols, smiles, formulas, parent, opers, weights, sascores,
            distToTarget, distToClosestDecoy, callerState, deliver);
        tbb::parallel_for(tbb::blocked_range<int>(0, morphAttempts),
            returnResults, tbb::auto_partitioner(), tbbCtx);
    }
    
    for (int i = 0; i < morphAttempts; ++i) {
        delete newMols[i];
    }
    delete[] newMols;
    delete[] opers;
    delete[] smiles;
    delete[] formulas;
    delete[] weights;
    delete[] sascores;
    delete[] distToTarget;
    delete[] distToClosestDecoy;

    delete mol;
    delete targetMol;
    delete targetFp;

    for (int i = 0; i < decoysFp.size(); ++i) {
        delete decoysFp[i];
    }

    for (int i = 0; i < strategies.size(); ++i) {
        delete strategies[i];
    }
    
}
