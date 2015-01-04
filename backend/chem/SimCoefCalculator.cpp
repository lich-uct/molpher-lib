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

#include <vector>

#include "chem/simCoefStrategy/SimCoefStrategy.h"
#include "chem/fingerprintStrategy/FingerprintStrategy.h"

// include strategies ..
#include "chem/fingerprintStrategy/AtomPairsFngpr.hpp"
#include "chem/fingerprintStrategy/MorganFngpr.hpp"
#include "chem/fingerprintStrategy/TopolLayeredFngpr1.hpp"
#include "chem/fingerprintStrategy/TopolLayeredFngpr2.hpp"
#include "chem/fingerprintStrategy/TopolSingleFngpr.hpp"
#include "chem/fingerprintStrategy/TopolTorsFngpr.hpp"
#include "chem/SimCoefCalculator.hpp"

// TODO: merge into one?
// include similarities 
#include "chem/simCoefStrategy/AllBitSimCoef.hpp"
#include "chem/simCoefStrategy/AsymmetricSimCoef.hpp"
#include "chem/simCoefStrategy/BraunBlanquetSimCoef.hpp"
#include "chem/simCoefStrategy/CosineSimCoef.hpp"
#include "chem/simCoefStrategy/DiceSimCoef.hpp"
#include "chem/simCoefStrategy/KulczynskiSimCoef.hpp"
#include "chem/simCoefStrategy/McConnaugheySimCoef.hpp"
#include "chem/simCoefStrategy/OnBitSimCoef.hpp"
#include "chem/simCoefStrategy/RusselSimCoef.hpp"
#include "chem/simCoefStrategy/SokalSimCoef.hpp"
#include "chem/simCoefStrategy/TanimotoSimCoef.hpp"
#include "chem/simCoefStrategy/TverskySimCoef.hpp"

// TODO: Remove use factory insted
SimCoefCalculator::SimCoefCalculator(
    SimCoeffSelector sc,
    FingerprintSelector fp,
    RDKit::ROMol *source,
    RDKit::ROMol *target
    )
{
    if (fp <= MAX_STANDARD_FP) {
        mExtended = false;
    } else {
        if (source == NULL || target == NULL) {
            mExtended = false;
        } else {
            mExtended = true;
        }
        fp = static_cast<FingerprintSelector>(fp - (MAX_STANDARD_FP + 1));
    }

    switch (sc) {
    case SC_TANIMOTO:
        mScStrategy = new TanimotoSimCoef();
        break;
    case SC_COSINE:
        mScStrategy = new CosineSimCoef();
        break;
    case SC_KULCZYNSKI:
        mScStrategy = new KulczynskiSimCoef();
        break;
    case SC_DICE:
        mScStrategy = new DiceSimCoef();
        break;
    case SC_TVERSKY_SUPERSTRUCTURE:
        /* Setting the weighting of prototype (in our case either target or
           neighborhood center) features to 90% (a=0.9) and variant (i.e. morph)
           features to 10% (b=0.1) means that mainly the prototype features
           are important, i.e., this produces a "superstucture-likeness" measure.
           In this case, a Tversky similarity value of 1.0 means that almost all
           prototype features are represented in the variant, 0.0 that almost
           none are. */
        mScStrategy = new TverskySimCoef(0.9, 0.1);
        break;
    case SC_TVERSKY_SUBSTRUCTURE:
        /* Conversely, setting the weights to 10% prototype (a=0.1) and
           90% variant (b=0.9) produces a "substructure-likeness" measure,
           where variants almost completely embedded into prototype have
           have values near 1.0. */
        mScStrategy = new TverskySimCoef(0.1, 0.9);
        break;
    case SC_SOKAL:
        mScStrategy = new SokalSimCoef();
        break;
    case SC_MC_CONNAUGHEY:
        mScStrategy = new McConnaugheySimCoef();
        break;
    case SC_ASYMMETRIC:
        mScStrategy = new AsymmetricSimCoef();
        break;
    case SC_BRAUN_BLANQUET:
        mScStrategy = new BraunBlanquetSimCoef();
        break;
    case SC_RUSSEL:
        mScStrategy = new RusselSimCoef();
        break;
    case SC_ON_BIT:
        mScStrategy = new OnBitSimCoef();
        break;
    case SC_ALL_BIT:
        mScStrategy = new AllBitSimCoef();
        break;
    default:
        mScStrategy = new TanimotoSimCoef();
        break;
    }

    switch (fp) {
    case FP_TOPOLOGICAL:
        mFpStrategy = new TopolSingleFngpr();
        break;
    case FP_TOPOLOGICAL_LAYERED_1:
        mFpStrategy = new TopolLayeredFngpr1();
        break;
    case FP_TOPOLOGICAL_LAYERED_2:
        mFpStrategy = new TopolLayeredFngpr2();
        break;
    case FP_ATOM_PAIRS:
        mFpStrategy = new AtomPairsFngpr();
        break;
    case FP_TOPOLOGICAL_TORSION:
        mFpStrategy = new TopolTorsFngpr();
        break;
    case FP_MORGAN:
        mFpStrategy = new MorganFngpr();
        break;
    default:
        mFpStrategy = new MorganFngpr();
        break;
    }

    if (!mExtended) {
        return;
    }

    // extended fingerprint
    std::set<AtomicNum> atomSet;

    AtomicNum atomicNum;
    for (int i = 0; i < source->getNumAtoms(); ++i) {
        atomicNum = source->getAtomWithIdx(i)->getAtomicNum();
        atomSet.insert(atomicNum);
    }
    for (int i = 0; i < target->getNumAtoms(); ++i) {
        atomicNum = target->getAtomWithIdx(i)->getAtomicNum();
        atomSet.insert(atomicNum);
    }

    int i = 0;
    for (std::set<AtomicNum>::iterator it = atomSet.begin();
            it != atomSet.end(); it++) {
        mAtomTypesToIdx[*it] = i;
        ++i;
    }
}

SimCoefCalculator::~SimCoefCalculator()
{
    delete mScStrategy;
    delete mFpStrategy;
}

double SimCoefCalculator::GetSimCoef(Fingerprint *fp1, Fingerprint *fp2)
{
    return mScStrategy->GetSimCoef(fp1, fp2);
}

double SimCoefCalculator::GetSimCoef(RDKit::ROMol *mol1, RDKit::ROMol *mol2)
{
    double result;
    Fingerprint *fp1 = GetFingerprint(mol1);
    Fingerprint *fp2 = GetFingerprint(mol2);
    result = mScStrategy->GetSimCoef(fp1, fp2);
    delete fp1;
    delete fp2;

    return result;
}

double SimCoefCalculator::GetSimCoef(Fingerprint *fp1, RDKit::ROMol *mol2)
{
    double result;
    Fingerprint *fp2 = GetFingerprint(mol2);
    result = mScStrategy->GetSimCoef(fp1, fp2);
    delete fp2;

    return result;
}

double SimCoefCalculator::ConvertToDistance(double coef) const
{
    return mScStrategy->ConvertToDistance(coef);
}

Fingerprint *SimCoefCalculator::GetFingerprint(RDKit::ROMol *mol)
{
    Fingerprint *fp = mFpStrategy->GetFingerprint(mol);
    if (mExtended) {
        Fingerprint *newFp = Extend(mol, fp);
        delete fp;
        fp = newFp;
    }
    return fp;
}

Fingerprint *SimCoefCalculator::Extend(RDKit::ROMol *mol, Fingerprint *fp)
{
    unsigned int maxBondOrder = RDKit::Bond::ONEANDAHALF;
    unsigned int atomTypes = mAtomTypesToIdx.size();
    unsigned int varSize = sizeof(extFpPart) * 8;

    int extraLength = (atomTypes + maxBondOrder + atomTypes * atomTypes) * varSize;
    int length = fp->getNumBits() + extraLength;

    Fingerprint *newFp = new Fingerprint(length);

    for (unsigned int i = 0; i < fp->getNumBits(); ++i) {
        if ((*fp)[i]) {
            newFp->setBit(i);
        }
    }

    std::vector<std::vector<extFpPart> > connectionTab;
	std::vector<extFpPart> cntBO;
	std::vector<extFpPart> cntAtomType;

    for (unsigned int i = 0; i < maxBondOrder; ++i) {
        cntBO.push_back(0);
    }

    for (unsigned int i = 0; i < atomTypes; ++i) {
		cntAtomType.push_back(0);

		std::vector<extFpPart> aux;
		for (unsigned int j = 0; j < atomTypes; ++j) {
            aux.push_back(0);
        }
		connectionTab.push_back(aux);
	}

	// atom types stats
	RDKit::Atom *atom;
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        atom = mol->getAtomWithIdx(i);
        ++cntAtomType[mAtomTypesToIdx[atom->getAtomicNum()]];
    }

	// bond stats
	RDKit::Bond *bond;
    for (unsigned int i = 0; i < mol->getNumBonds(); ++i) {
        bond = mol->getBondWithIdx(i);
		AtomicNum typeBegin = bond->getBeginAtom()->getAtomicNum();
        AtomicNum typeEnd = bond->getEndAtom()->getAtomicNum();
		int bo = static_cast<int>(bond->getBondType());

		connectionTab[mAtomTypesToIdx[typeBegin]][mAtomTypesToIdx[typeEnd]] += bo;
		connectionTab[mAtomTypesToIdx[typeEnd]][mAtomTypesToIdx[typeBegin]] += bo;
		if (bo < maxBondOrder) {
            ++cntBO[bo];
        }
	}
	for (unsigned int i = 0; i < atomTypes; ++i) {
        connectionTab[i][i] /= 2;
    }

	// insert the stats into the FP
    unsigned int idx = fp->getNumBits();
    extFpPart mask;

    for (unsigned int i = 0; i < atomTypes; ++i) {
        mask = 1;
        for (unsigned int sh = 0; sh < varSize; ++sh) {
            if (cntAtomType[i] & (mask << sh)) {
                newFp->setBit(idx);
            }
            ++idx;
        }
    }

    for (unsigned int i = 0; i < maxBondOrder; ++i) {
        mask = 1;
        for (unsigned int sh = 0; sh < varSize; ++sh) {
            if (cntBO[i] & (mask << sh)) {
                newFp->setBit(idx);
            }
            ++idx;
        }
    }

    for (unsigned int i = 0; i < atomTypes; ++i) {
        for (unsigned int j = 0; j < atomTypes; ++j) {
            mask = 1;
            for (unsigned int sh = 0; sh < varSize; ++sh) {
                if (connectionTab[i][j] & (mask << sh)) {
                    newFp->setBit(idx);
                }
                ++idx;
            }
        }
    }

    return newFp;
}
