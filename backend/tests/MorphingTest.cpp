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

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <tbb/task.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

#include "inout.h"
#include "auxiliary/SynchRand.h"
#include "chem/fingerprintStrategy/FingerprintStrategy.h"
#include "chem/ChemicalAuxiliary.h"
#include "chem/simCoefStrategy/SimCoefStrategy.h"
#include "chem/morphing/Morphing.hpp"
#include "MorphingTest.h"
#include "extensions/SAScore.h"
#include "chem/SimCoefCalculator.hpp"

using namespace std;

// testing the modification of the formal charge
void TestSetFormalCharge(RDKit::RWMol &mol)
{
    cout << "Test formal charge" << endl;

    int i = 0;
    RDKit::Atom *atom;
    for (RDKit::ROMol::AtomIterator iter = mol.beginAtoms();
            iter != mol.endAtoms(); iter++, i++) {
        atom = *iter;
        cout << "the " << i << "th atom is " << atom->getSymbol() <<
            " the formal charge is" << atom->getFormalCharge() << endl;
    }

    SetFormalCharge(1, mol);

    i = 0;
    for (RDKit::ROMol::AtomIterator iter = mol.beginAtoms();
            iter != mol.endAtoms(); iter++, i++) {
        atom = *iter;
        cout << "the " << i << "th atom is " << atom->getSymbol() <<
            " the formal charge is" << atom->getFormalCharge() << endl;
    }
}

void TestSerialization(RDKit::ROMol &mol)
{
    cout << "serialization string:" << endl;
    string pickle;
    RDKit::MolPickler::pickleMol(mol, pickle);

    cout << pickle << endl;
    cout << "size of the pickle is " << pickle.length() << endl;
    cout << "size of the molecule is " << sizeof(mol) << endl;
}

// writes important informations about i-th atom of the molecule
void WriteAtomInfo(RDKit::ROMol *mol, int i)
{
    RDKit::Atom *atom = mol->getAtomWithIdx(i);
    atom->calcExplicitValence();
    atom->calcImplicitValence();
    RDKit::PeriodicTable *table = RDKit::PeriodicTable::getTable();
    cout << "atom " << i << endl;
    cout << "   has atomic number " << atom->getAtomicNum() << endl;
    cout << "   symbol " << atom->getSymbol() << endl;
    cout << "   degree " << atom->getDegree() << endl;
    cout << "   total number Hs " << atom->getTotalNumHs() << endl;
}

void WriteAtomInfoFromMolecule(RDKit::ROMol *mol)
{
    for (int i = 0; i < mol->getNumAtoms(); ++i) {
        WriteAtomInfo(mol, i);
    }
}

double TestSimCoef(RDKit::ROMol *mol1, RDKit::ROMol *mol2)
{
    SimCoefCalculator sim(SC_TANIMOTO, FP_TOPOLOGICAL);
    return sim.GetSimCoef(mol1, mol2);
}

RDKit::RWMol *CreatAlcohol()
{
    RDKit::RWMol *result =  new RDKit::RWMol();

    // create new atoms
    RDKit::Atom *at0 = new RDKit::Atom(6); // C
    RDKit::Atom *at1 = new RDKit::Atom(8); // O

    // add atoms to molecule
    result->addAtom(at0); // idx in mol 0 - C
    result->addAtom(at0); // idx in mol 1 - C
    result->addAtom(at1); // idx in mol 2 - O

    // delete atoms
    delete at0;
    delete at1;

    // add bound to atom
    result->addBond(result->getAtomWithIdx(0), result->getAtomWithIdx(1),
        RDKit::Bond::SINGLE); // C1-C2
    result->addBond(result->getAtomWithIdx(1), result->getAtomWithIdx(2),
        RDKit::Bond::SINGLE); // C2-O

    return result;
}

RDKit::RWMol *CreatAlcoholWitH()
{
    RDKit::RWMol *result =  new RDKit::RWMol();

    // create new atoms
    RDKit::Atom *at0 = new RDKit::Atom(6); // C
    RDKit::Atom *at1 = new RDKit::Atom(8); // O
    RDKit::Atom *at2 = new RDKit::Atom(1); // H

    // add atoms to molecule
    result->addAtom(at0); // idx in mol 0 - C
    result->addAtom(at0); // idx in mol 1 - C
    result->addAtom(at1); // idx in mol 2 - O
    result->addAtom(at2); // idx in mol 3 - H
    result->addAtom(at2); // idx in mol 4 - H
    result->addAtom(at2); // idx in mol 5 - H
    result->addAtom(at2); // idx in mol 6 - H
    result->addAtom(at2); // idx in mol 7 - H
    result->addAtom(at2); // idx in mol 8 - H

    // delete atoms
    delete at0;
    delete at1;
    delete at2;

    // add bound to atom
    result->addBond(result->getAtomWithIdx(0), result->getAtomWithIdx(1),
        RDKit::Bond::SINGLE); // C1-C2
    result->addBond(result->getAtomWithIdx(1), result->getAtomWithIdx(2),
        RDKit::Bond::SINGLE); // C2-O
    result->addBond(result->getAtomWithIdx(0), result->getAtomWithIdx(3)); // C1-H1
    result->addBond(result->getAtomWithIdx(0), result->getAtomWithIdx(4)); // C1-H2
    result->addBond(result->getAtomWithIdx(0), result->getAtomWithIdx(5)); // C1-H3
    result->addBond(result->getAtomWithIdx(1), result->getAtomWithIdx(6)); // C2-H4
    result->addBond(result->getAtomWithIdx(1), result->getAtomWithIdx(7)); // C2-H5
    result->addBond(result->getAtomWithIdx(2), result->getAtomWithIdx(8)); // O-H6

    return result;
}

void TestAlcohol()
{
    RDKit::RWMol *mol = CreatAlcohol();
    RDKit::RWMol *sMol = CreatAlcohol();
    RDKit::MolOps::sanitizeMol(*sMol);

    cout << "alcohol info: " << endl;
    WriteAtomInfoFromMolecule(mol);

    cout << "alcohol info with H: " << endl;
    WriteAtomInfoFromMolecule(sMol);

    cout << "alcohol diff: ";
    cout << TestSimCoef(mol, sMol) << endl;
}

void TestBondBetween()
{
    RDKit::RWMol *mol =  new RDKit::RWMol();

    // create new atoms
    RDKit::Atom *at0 = new RDKit::Atom(6); // C
    RDKit::Atom *at1 = new RDKit::Atom(8); // O

    // add atoms to molecule
    mol->addAtom(at0); // idx in mol 0 - C
    mol->addAtom(at0); // idx in mol 1 - C
    mol->addAtom(at1); // idx in mol 2 - O

    // delete atoms
    delete at0;
    delete at1;

    // add bound to atom
    mol->addBond(mol->getAtomWithIdx(0), mol->getAtomWithIdx(1),
        RDKit::Bond::SINGLE); // C1-C2
    mol->addBond(mol->getAtomWithIdx(1), mol->getAtomWithIdx(2),
        RDKit::Bond::SINGLE); // C2-O

    RDKit::Bond *bond;

    if (bond = mol->getBondBetweenAtoms(1, 2)) {
        cout << "there is bond between 1-2" << endl;
        cout << "    " <<  bond->getOtherAtom(bond->getBeginAtom())->getIdx() << endl;
        cout << "    " <<  bond->getOtherAtom(bond->getEndAtom())->getIdx() << endl;
    }

    if (bond = mol->getBondBetweenAtoms(2, 1)) {
        cout << "there is bond between 2-1" << endl;
        cout << "    " << bond->getOtherAtom(bond->getBeginAtom())->getIdx() << endl;
        cout << "    " << bond->getOtherAtom(bond->getEndAtom())->getIdx() << endl;
    }

    if (bond = mol->getBondBetweenAtoms(0, 2)) {
        cout << "there is bond between 0-2" << endl;
    } else {
        cout << "there is no bond between 0-2" << endl;
    }

    delete mol;
}

void SerializationBenchmark(RDKit::ROMol &mol)
{
    int iterations = 10000;
    string str;
    clock_t start;
    clock_t finish;

    cout << "Pickle benchmark:" << endl;
    RDKit::MolPickler::pickleMol(mol, str);
    cout << "Pickle size [Bytes] = " << str.length() << endl;
    start = clock();
    for (int i = 0; i < iterations; ++i) {
        string pickle;
        RDKit::MolPickler::pickleMol(mol, pickle);
        RDKit::ROMol *molecule = new RDKit::ROMol();
        RDKit::MolPickler::molFromPickle(pickle, molecule);
        delete molecule;
    }
    finish = clock();
    cout << "Pickle time [msec] = " << finish - start << endl;

    str.clear();
    cout << endl;

    cout << "Smiles benchmark:" << endl;
    str = RDKit::MolToSmiles(mol);
    cout << "Smiles size [Bytes] = " << str.length() << endl;
    start = clock();
    for (int i = 0; i < iterations; ++i) {
        string smile = RDKit::MolToSmiles(mol);
        RDKit::ROMol *molecule = RDKit::SmilesToMol(smile);
        delete molecule;
    }
    finish = clock();
    cout << "Smiles time [msec] = " << finish - start << endl;
}

void FPandSCBenchmark(RDKit::ROMol *mol1, RDKit::ROMol *mol2)
{
    int iterations = 10000;
    clock_t start;
    clock_t finish;

    for (int i = 0; i <= FP_MORGAN; ++i) {
        FingerprintSelector fp = static_cast<FingerprintSelector>(i);
        SimCoefCalculator sCC(SC_TANIMOTO, fp);
        Fingerprint **fps = new Fingerprint *[iterations];

        start = clock();
        for (int j = 0; j < iterations; ++j) {
            fps[j] = sCC.GetFingerprint(mol1);
        }
        finish = clock();
        cout << "FPSelector " << fp << " time [msec] = " <<
            finish - start << endl;

        for (int j = 0; j < iterations; ++j) {
            delete fps[j];
        }
        delete[] fps;
    }
    cout << endl;

    for (int i = 0; i <= SC_ALL_BIT; ++i) {
        SimCoeffSelector sc = static_cast<SimCoeffSelector>(i);
        Fingerprint *fp1 = GetFingerprint(mol1, FP_MORGAN);
        Fingerprint *fp2 = GetFingerprint(mol2, FP_MORGAN);
        SimCoefCalculator sCC(sc, FP_MORGAN);
        start = clock();
        for (int j = 0; j < iterations; ++j) {
            sCC.GetSimCoef(fp1, fp2);
        }
        finish = clock();
        cout << "SCSelector " << sc << " time [msec] = " <<
            finish - start << endl;
        delete fp1;
        delete fp2;
    }
}

void MolToMolBlockBenchmark(RDKit::ROMol *mol)
{
    int iterations = 10000;
    clock_t start;
    clock_t finish;

    string molBlock;
    RDKit::ROMol **newMol = new RDKit::ROMol *[iterations];

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        molBlock = RDKit::MolToMolBlock(*mol);
        newMol[i] = RDKit::MolBlockToMol(molBlock);
    }
    finish = clock();
    cout << "mol => molBlock => mol time [msec] = " <<
        finish - start << endl;

    for (int i = 0; i < iterations; ++i) {
        delete newMol[i];
    }
    delete[] newMol;
}

// benchmark generation of molecule formula against SMILE generation
void BenchmarkMolBlockVsSmiles(RDKit::ROMol *mol)
{
    int iterations = 10000;
    clock_t start;
    clock_t finish;

    string strMol;
    RDKit::ROMol **newMol = new RDKit::ROMol *[iterations];

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        strMol = RDKit::MolToMolBlock(*mol);
        newMol[i] = RDKit::MolBlockToMol(strMol);
    }
    finish = clock();
    cout << "mol => molBlock => mol time [msec] = " <<
        finish - start << endl;
    cout << "length = " << strMol.length() << endl;

    for (int i = 0; i < iterations; ++i) {
        delete newMol[i];
    }

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        strMol = RDKit::MolToSmiles(*mol);
        newMol[i] = RDKit::SmilesToMol(strMol);
        if (i == 0) cout << strMol.length() << endl;
    }
    finish = clock();
    cout << "mol => smiles => mol time [msec] = " <<
        finish - start << endl;
    cout << "length = " << strMol.length() << endl;


    for (int i = 0; i < iterations; ++i) {
        delete newMol[i];
    }

    delete[] newMol;
}

void BenchmarkSmilesVsMolFormula(RDKit::ROMol *mol)
{
    int iterations = 10000;
    clock_t start;
    clock_t finish;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        RDKit::MolToSmiles(*mol);
    }
    finish = clock();
    cout << "mol => Smiles time [msec] = " << finish - start << endl;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        RDKit::Descriptors::calcMolFormula(*mol);
    }
    finish = clock();
    cout << "mol => mol formula time [msec] = " << finish - start << endl;
}

void PrintNeighbors(RDKit::Atom &atom, std::ostream &st)
{
    if (atom.getIsAromatic()) {
        st << "the atom is aromatic" << std::endl;
    }
    std::string nAtName;
    AtomIdx nAtomIdx;
    RDKit::Bond::BondType bt;
    RDKit::ROMol &mol = atom.getOwningMol();
    RDKit::ROMol::ADJ_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomNeighbors(&atom);
    while (beg != end) {
        nAtomIdx = (mol)[*beg].get()->getIdx();
        nAtName = (mol)[*beg].get()->getSymbol();
        bt = mol.getBondBetweenAtoms(atom.getIdx(), nAtomIdx)->getBondType();
        st << atom.getIdx() << "(" << atom.getSymbol() << ")" << "-"
            << nAtomIdx << "(" << nAtName << ")" << ", bond order " << bt << std::endl;
        beg++;
    }
    st << std::endl;
}

// tests the valence of the atoms in the molecule, if there is a problem
// writes the bonds of the atom
void TestMol(RDKit::RWMol &mol, string err, ostream &st)
{
    RDKit::Atom *atom;
    for (int i = 0; i < mol.getNumAtoms(); ++i)
    {
        atom = mol.getAtomWithIdx(i);
        try {
            atom->calcExplicitValence();
            atom->calcImplicitValence();
        } catch(exception &e) {
            st << "error occured at " << err << endl;
            PrintNeighbors(*atom, st);
        }
    }
}

// the same function as RDKit::MolOps::sanitizeMol, but there is molecule checking
void MySanitize(RDKit::RWMol &mol, ostream &st)
{
    try {
        mol.clearComputedProps();
    } catch (exception &e) {
        TestMol(mol, "clearComputedProps", st);
    }

    try {
        RDKit::MolOps::cleanUp(mol);
    } catch (exception &e) {
        TestMol(mol, "cleanUp", st);
    }

    try {
        mol.updatePropertyCache();
    } catch (exception &e) {
        TestMol(mol, "updatePropertyCache", st);
    }

    try {
        RDKit::MolOps::Kekulize(mol);
    } catch (exception &e) {
        TestMol(mol, "Kekulize", st);
    }

    try {
        RDKit::MolOps::setAromaticity(mol);
    } catch (exception &e) {
        TestMol(mol, "setAromaticity", st);
    }

    try {
        RDKit::MolOps::setConjugation(mol);
    } catch (exception &e) {
        TestMol(mol, "setConjugation", st);
    }

    try {
        RDKit::MolOps::setHybridization(mol);
    } catch (exception &e) {
        TestMol(mol, "setHybridization", st);
    }

    try {
        RDKit::MolOps::cleanupChirality(mol);
    } catch (exception &e) {
        TestMol(mol, "cleanupChirality", st);
    }

    try {
        RDKit::MolOps::adjustHs(mol);
    } catch (exception &e) {
        TestMol(mol, "adjustHs", st);
    }
}

void TestInterlay(RDKit::RWMol mol, AtomIdx at1, AtomIdx at2, AtomicNum newAt, ostream &st)
{
    st << "molecule before modifying" << endl;
    st << RDKit::MolToMolBlock(mol) << endl << endl;
    st << "bond to interly: " << at1 + 1 << "-" << at2 + 1 << endl << endl;

    RDKit::Atom atom(newAt);
    double mass =
        RDKit::PeriodicTable::getTable()->getMostCommonIsotopeMass(atom.getAtomicNum());
    atom.setMass(mass);

    RDKit::Bond::BondType bt = mol.getBondBetweenAtoms(0, 15)->getBondType();
    int newAtom = mol.addAtom(&atom);
    st << "the new atom " << newAtom + 1 << endl;
    mol.addBond(0, newAtom, bt);
    mol.addBond(15, newAtom, bt);
    mol.removeBond(0, 15);

//    try {
//        string smiles = RDKit::MolToSmiles(mol);
//        RDKit::SmilesToMol(smiles);
//    } catch (exception &e) {
//        cout << "some error at mol => smiles => mol";
//    }

//    MySanitize(mol, st);

    try {
        st << RDKit::MolToMolBlock(mol) << endl << endl;
//        string smiles = RDKit::MolToSmiles(mol);
//        RDKit::ROMol newMol = *(RDKit::SmilesToMol(smiles));
//        st << RDKit::MolToMolBlock(newMol);
    } catch (exception &e) {
        cout << "some error after sanitize at mol => smiles => mol";
    }
}

void TestWriteToSDF(std::vector<RDKit::RWMol *> &mols)
{
    RDKit::SDWriter writer("sdfile.sdf");
    for (int i = 0; i < mols.size(); ++i) {
        writer.write(*(mols[i]));
        writer.flush();
    }
    writer.close();
}

// replaces aromatic bonds with single and double bonds
void RecalculateAromaticRings(RDKit::RWMol &mol)
{
    RDKit::VECT_INT_VECT ringVect = mol.getRingInfo()->atomRings();

    RDKit::Bond::BondType bt;
    for (int i = 0; i < ringVect.size(); ++i) {
        // if the ring is aromatic
        if (mol.getAtomWithIdx(ringVect[i][0])->getIsAromatic()) {
            // decides whether aromatic ring starts with single or double bond
            int rand = SynchRand::GetRandomNumber(0, 1);
            for (int j = 0; j < ringVect[i].size(); ++j) {
                if (j % 2 == rand) {
                    bt = RDKit::Bond::SINGLE;
                } else {
                    bt = RDKit::Bond::DOUBLE;
                }

                AtomIdx begin = ringVect[i][j];
                AtomIdx end = ringVect[i][(j + 1) % ringVect[i].size()];

                mol.getAtomWithIdx(begin)->setIsAromatic(false);
                mol.getAtomWithIdx(end)->setIsAromatic(false);
                mol.getBondBetweenAtoms(begin, end)->setBondType(bt);
                mol.getBondBetweenAtoms(begin, end)->setIsAromatic(false);
            }
        }
    }

    RDKit::Atom *atom;
    for (int i = 0; i < mol.getNumAtoms(); ++i) {
        atom = mol.getAtomWithIdx(i);
        atom->calcExplicitValence();
        atom->calcImplicitValence();
    }
}

void TestRingReset(RDKit::RWMol mol)
{
    RecalculateAromaticRings(mol);

    mol.debugMol(cout);
    {
        RDKit::VECT_INT_VECT ringVect = mol.getRingInfo()->atomRings();

        RDKit::Bond::BondType bt;
        for (int i = 0; i < ringVect.size(); ++i) {
            cout << "ring " << i+1 << ":" << endl;
            for (int j = 0; j < ringVect[i].size(); ++j) {
                AtomIdx begin = ringVect[i][j];
                AtomIdx end = ringVect[i][(j + 1) % ringVect[i].size()];
                bt = mol.getBondBetweenAtoms(begin, end)->getBondType();
                cout << begin << "-" << end << ", bond type: " << bt << endl;
            }
        }
    }

    mol.getRingInfo()->reset();

    {
        RDKit::VECT_INT_VECT ringVect = mol.getRingInfo()->atomRings();

        RDKit::Bond::BondType bt;
        for (int i = 0; i < ringVect.size(); ++i) {
            cout << "ring " << i+1 << ":" << endl;
            for (int j = 0; j < ringVect[i].size(); ++j) {
                AtomIdx begin = ringVect[i][j];
                AtomIdx end = ringVect[i][(j + 1) % ringVect[i].size()];
                bt = mol.getBondBetweenAtoms(begin, end)->getBondType();
                cout << begin << "-" << end << ", bond type: " << bt << endl;
            }
        }
    }

    RDKit::MolOps::sanitizeMol(mol);

    {
        RDKit::VECT_INT_VECT ringVect = mol.getRingInfo()->atomRings();

        RDKit::Bond::BondType bt;
        for (int i = 0; i < ringVect.size(); ++i) {
            cout << "ring " << i+1 << ":" << endl;
            for (int j = 0; j < ringVect[i].size(); ++j) {
                AtomIdx begin = ringVect[i][j];
                AtomIdx end = ringVect[i][(j + 1) % ringVect[i].size()];
                bt = mol.getBondBetweenAtoms(begin, end)->getBondType();
                cout << begin << "-" << end << ", bond type: " << bt << endl;
            }
        }
    }
}

void DummyDeliver(MolpherMolecule *mol, void *state)
{
    // no-op
}

void TestMorphing(RDKit::ROMol *source, RDKit::ROMol *target)
{
    string sSmile = RDKit::MolToSmiles(*source);
    string tSmile = RDKit::MolToSmiles(*target);

    MolpherMolecule sMol(sSmile);
    MolpherMolecule tMol(tSmile);

    std::vector<ChemOperSelector> operators;
    operators.push_back(OP_ADD_ATOM);
    operators.push_back(OP_REMOVE_ATOM);
    operators.push_back(OP_ADD_BOND);
    operators.push_back(OP_REMOVE_BOND);
    operators.push_back(OP_MUTATE_ATOM);
    operators.push_back(OP_INTERLAY_ATOM);
    operators.push_back(OP_BOND_REROUTE);
    operators.push_back(OP_BOND_CONTRACTION);

    vector<MolpherMolecule> decoys;
    tbb::task_group_context tbbCtx;
    GenerateMorphs(sMol, 5000, FP_MORGAN, SC_TANIMOTO,
        operators, tMol, decoys, tbbCtx, NULL, DummyDeliver);
}

void TestRemoveRing(RDKit::RWMol mol)
{
    RecalculateAromaticRings(mol);

    mol.debugMol(cout);

    RDKit::VECT_INT_VECT ringVect = mol.getRingInfo()->atomRings();

    RDKit::Bond::BondType bt;
    for (int i = 0; i < ringVect.size(); ++i) {
        cout << "ring " << i + 1 << ":" << endl;
        for (int j = 0; j < ringVect[i].size(); ++j) {
            AtomIdx begin = ringVect[i][j];
            AtomIdx end = ringVect[i][(j + 1) % ringVect[i].size()];
            bt = mol.getBondBetweenAtoms(begin, end)->getBondType();
            cout << begin << "-" << end << ", bond type: " << bt << endl;
        }
    }

    cout << endl << endl << "delete bonds" << endl;
    for (int i = 0; i < ringVect.size(); ++i) {
        for (int j = 0; j < ringVect[i].size(); ++j) {
            RDKit::RWMol newMol(mol);
            AtomIdx begin = ringVect[i][j];
            AtomIdx end = ringVect[i][(j + 1) % ringVect[i].size()];
            newMol.removeBond(begin, end);

            cout << begin << "-" << end << endl;
            PrintNeighbors(*(newMol.getAtomWithIdx(begin)), cout);
            PrintNeighbors(*(newMol.getAtomWithIdx(end)), cout);

            string smile = RDKit::MolToSmiles(newMol);
            RDKit::SmilesToMol(smile);

            MySanitize(newMol, cout);
        }
        cout << endl << endl;
    }
}

void TestTransform(RDKit::RWMol mol)
{
    fstream st;

    cout << RDKit::MolToSmiles(mol) << endl;

    st.open ("sdf/input.sdf", fstream::out);
    st << RDKit::MolToMolBlock(mol) << endl;
    st.close();

    RDKit::RWMol query = *RDKit::SmilesToMol("c1ccc(C#N)c(C(F)(F)F)c1");
    st.open ("sdf/query.sdf", fstream::out);
    st << RDKit::MolToMolBlock(query) << endl;
    st.close();

    RDKit::RWMol replacement = *RDKit::SmilesToMol("C=C");
    st.open ("sdf/replacement.sdf", fstream::out);
    st << RDKit::MolToMolBlock(replacement) << endl;
    st.close();

    vector<RDKit::ROMOL_SPTR> vect = RDKit::replaceSubstructs(mol, query, replacement);

    for (int i = 0; i < vect.size(); ++i) {
        stringstream ss;
        ss << "sdf/output" << i << ".sdf";
        st.open (ss.str().c_str(), fstream::out);
        RDKit::ROMol result = *vect[i];
        result.debugMol(cout);
        try {
            st << RDKit::MolToMolBlock(result) << endl;
        } catch (ValueErrorException &e) {
            cout << e.message() << endl;
        }
        st.close();
    }
}

void SimpleTest(RDKit::RWMol mol, BondIdx bidx)
{
    RDKit::MolOps::Kekulize(mol);
    mol.getRingInfo()->reset();

    AtomIdx begin = mol.getBondWithIdx(bidx)->getBeginAtomIdx();
    AtomIdx end = mol.getBondWithIdx(bidx)->getEndAtomIdx();

    mol.removeBond(begin, end);

    try {
        RDKit::RWMol newMol;
        CopyMol(mol, newMol);

        try {
            RDKit::MolOps::sanitizeMol(newMol);
            fstream st;
            stringstream ss;
            ss << "sdf/newMol" << bidx << ".sdf";
            st.open (ss.str().c_str(), fstream::out);
            st << RDKit::MolToMolBlock(newMol) << endl;
            st.close();
            cout << bidx << " ok" << endl;
        } catch (ValueErrorException &e) {
            cout << bidx << " error" << endl;
            cout << e.message() << endl;
        } catch (exception &e) {
            cout << bidx << " error" << endl;
            cout << e.what() << endl;
        }
    } catch (ValueErrorException &e) {
        cout << bidx << " error" << endl;
        cout << e.message() << endl;
    } catch (exception &e) {
        cout << bidx << " error" << endl;
        cout << e.what() << endl;
    }

//    mol.getRingInfo()->initialize();
//    try {
//        RDKit::MolOps::sanitizeMol(mol);
//        mol.debugMol(cout);
//    } catch (ValueErrorException &e) {
//        cout << e.message() << endl;
//    } catch (exception &e) {
//        cout << e.what() << endl;
//    }
}

void SimpleTest(RDKit::RWMol mol)
{
    RDKit::Bond *bond;
    for (int i = 0; i < mol.getNumBonds(); ++i) {
        bond = mol.getBondWithIdx(i);
        if (RDKit::queryIsBondInRing(bond)) {
            SimpleTest(mol, i);
        }
    }
}

void TestInterlay(RDKit::RWMol mol, AtomIdx at1, AtomIdx at2, AtomicNum newAt)
{
    RDKit::MolOps::Kekulize(mol);
    mol.getRingInfo()->reset();

    RDKit::Atom atom(newAt);
    double mass =
        RDKit::PeriodicTable::getTable()->getMostCommonIsotopeMass(atom.getAtomicNum());
    atom.setMass(mass);

    RDKit::Bond::BondType bt = mol.getBondBetweenAtoms(0, 15)->getBondType();
    int newAtom = mol.addAtom(&atom);
    mol.addBond(0, newAtom, bt);
    mol.addBond(15, newAtom, bt);
    mol.removeBond(0, 15);

    try {
        RDKit::RWMol newMol;
        CopyMol(mol, newMol);

        try {
            RDKit::MolOps::sanitizeMol(newMol);
            fstream st;
            stringstream ss;
            ss << "sdf/interlay" << at1 << "-" << at2 << ".sdf";
            st.open (ss.str().c_str(), fstream::out);
            st << RDKit::MolToMolBlock(newMol) << endl;
            st.close();
            cout << at1 << "-" << at2 << " ok" << endl;
        } catch (ValueErrorException &e) {
            cout << at1 << "-" << at2 << " error" << endl;
            cout << e.message() << endl;
        } catch (exception &e) {
            cout << at1 << "-" << at2 << " error" << endl;
            cout << e.what() << endl;
        }
    } catch (ValueErrorException &e) {
        cout << at1 << "-" << at2 << " error" << endl;
        cout << e.message() << endl;
    } catch (exception &e) {
        cout << at1 << "-" << at2 << " error" << endl;
        cout << e.what() << endl;
    }

//    mol.getRingInfo()->initialize();
//    try {
//        RDKit::MolOps::sanitizeMol(mol);
//        mol.debugMol(cout);
//    } catch (ValueErrorException &e) {
//        cout << e.message() << endl;
//    } catch (exception &e) {
//        cout << e.what() << endl;
//    }
}

void DictBenchmark()
{
    vector<string> data;
    data.push_back("string1");
    data.push_back("string2");
    data.push_back("string3");
    data.push_back("string4");
    data.push_back("string5");
    data.push_back("string6");
    data.push_back("string7");
    data.push_back("string8");
    data.push_back("string9");
    data.push_back("string10");
    data.push_back("string11");
    data.push_back("string12");
    data.push_back("string13");
    data.push_back("string14");
    data.push_back("string15");
    data.push_back("string16");
    data.push_back("string17");
    data.push_back("string18");
    data.push_back("string19");
    data.push_back("string20");

    int iterations = 100000;
    clock_t start;
    clock_t finish;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        vector<string> copy;
        copy = data;
    }
    finish = clock();
    cout << "Copy ctor = " << finish - start << " msec" << endl;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        stringstream txtStream(ios_base::out | ios_base::in);
        boost::archive::text_oarchive txtOutArchive(txtStream);
        txtOutArchive << data;
        string result = txtStream.str();
    }
    finish = clock();
    cout << "Serialization (text out) = " << finish - start << " msec" << endl;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        stringstream binStream(ios_base::out | ios_base::in | ios_base::binary);
        boost::archive::binary_oarchive binOutArchive(binStream);
        binOutArchive << data;
        string result = binStream.str();
    }
    finish = clock();
    cout << "Serialization (bin out) = " << finish - start << " msec" << endl;

    string stored;

    {
        stringstream txtStream(ios_base::out | ios_base::in);
        boost::archive::text_oarchive txtOutArchive(txtStream);
        txtOutArchive << data;
        stored = txtStream.str();
    }

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        vector<string> restored;
        stringstream txtStream(stored, ios_base::out | ios_base::in);
        boost::archive::text_iarchive txtInArchive(txtStream);
        txtInArchive >> restored;
    }
    finish = clock();
    cout << "Serialization (text in) = " << finish - start << " msec" << endl;

    {
        stringstream binStream(ios_base::out | ios_base::in | ios_base::binary);
        boost::archive::binary_oarchive binOutArchive(binStream);
        binOutArchive << data;
        stored = binStream.str();
    }

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        vector<string> restored;
        stringstream binStream(stored, ios_base::out | ios_base::in | ios_base::binary);
        boost::archive::binary_iarchive binInArchive(binStream);
        binInArchive >> restored;
    }
    finish = clock();
    cout << "Serialization (bin in) = " << finish - start << " msec" << endl;
}

void CopyBenchmark(RDKit::RWMol &mol)
{
    int iterations = 100000;
    clock_t start;
    clock_t finish;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        string molBlock = RDKit::MolToMolBlock(mol);
        RDKit::RWMol newMol = *RDKit::MolBlockToMol(molBlock);
    }
    finish = clock();
    cout << "mol => molBlock => mol time [msec] = " <<
        finish - start << endl;

    start = clock();
    for (int i = 0; i < iterations; ++i) {
        RDKit::RWMol newMol;
        CopyMol(mol, newMol);
    }
    finish = clock();
    cout << "CopyMol time [msec] = " <<
        finish - start << endl;
}

void TestSAScore(RDKit::ROMol *mol)
{
    SAScore::loadData();
    SAScore* score = SAScore::getInstance();
    std::cout << score->getScore(*mol);
    SAScore::destroyInstance();
}

void MorphingSandbox()
{
//    DictBenchmark();
//    return;

    string path = "TestFiles/CID_10635-CID_15951529.sdf";
    std::vector<RDKit::RWMol *> mols;
    ReadRWMolsFromSDF(path, mols);
//    CopyBenchmark(*mols[1]);
//    return;
//    TestTransform(*mols[1]);
//    SimpleTest(*mols[1]);
//    TestInterlay(*mols[1], 0, 15, 6);
//    return;
//    BenchmarkMolBlockVsSmiles(mols[1]);
//    return;

//    TestRingReset(*mols[1]);
//    return;

//    TestRemoveRing(*mols[1]);
//    return;
//
//    TestMorphing(mols[0], mols[1]);
    RDKit::RWMol * molecule = RDKit::SmilesToMol ("CNC1CC(=O)C(F)=CC1(C)N1C2CN(O)C(C)(C2)C(C)C1", 0, true, 0);
    TestSAScore(molecule);
    delete molecule;
//    TestMorphing(mols[1], mols[0]);
    return;

//    SerializationBenchmark(*mols[0]);
//    FPandSCBenchmark(mols[0], mols[1]);
//    MolToMolBlockBenchmark(mols[1]);

    BenchmarkMolBlockVsSmiles(mols[1]);
    BenchmarkSmilesVsMolFormula(mols[1]);
    return;

//    TestWriteToSDF(mols);
//    return;


    RDKit::RWMol mol = RDKit::RWMol(*mols[1]);

//    cout << RDKit::MolToMolBlock(mol) << endl;
//    RDKit::MolOps::sanitizeMol(mol);
//    cout << RDKit::MolToMolBlock(mol) << endl;
//    return;

    fstream st;
    st.open("errors.txt", fstream::out);

    bool closed = false;
    try {
        TestInterlay(RDKit::RWMol(*mols[1]), 0, 15, 6, st);
        closed = true;
        st.close();
    } catch (exception &e) {
        if (!closed) {
            st.close();
        }
    }

    return;
}


