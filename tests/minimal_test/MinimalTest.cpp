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

#include <stdexcept>
#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <random_seed.hpp>
#include <morphing/operators/AddAtom.hpp>
#include <morphing/Molpher.hpp>
#include <mol_helpers.hpp>
#include <morphing/operators/RemoveAtom.hpp>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <morphing/operators/AddBond.hpp>
#include <morphing/operators/RemoveBond.hpp>
#include <morphing/operators/MutateAtom.hpp>
#include <morphing/operators/InterlayAtom.hpp>
#include <morphing/operators/ContractBond.hpp>
#include <morphing/operators/RerouteBond.hpp>

#include "MinimalTest.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(MinimalTest);

MinimalTest::MinimalTest() :
test_dir("tests/test_files/")
{
    // no action
}

MinimalTest::~MinimalTest() {
    // no action
}

void MinimalTest::setUp() {
    load_data_from("src/python/molpher/SAScore.dat");
    set_random_seed(42);
}

void MinimalTest::tearDown() {
    // no action
}

void MinimalTest::testAtom() {
	MolpherAtom carbon = MolpherAtom("C");
	MolpherAtom oxygen = MolpherAtom("O");
	CPPUNIT_ASSERT_EQUAL(carbon.getFormalCharge(), 0);
	CPPUNIT_ASSERT_EQUAL(oxygen.getFormalCharge(), 0);
	oxygen.setFormalCharge(-1);
	carbon.setFormalCharge(+1);
	CPPUNIT_ASSERT_EQUAL(carbon.getFormalCharge(), +1);
	CPPUNIT_ASSERT_EQUAL(oxygen.getFormalCharge(), -1);

	// locking features
	oxygen.setLockingMask(MolpherAtom::NO_ADDITION | MolpherAtom::NO_MUTATION);
	CPPUNIT_ASSERT(oxygen.isLocked());
	CPPUNIT_ASSERT(MolpherAtom::NO_ADDITION & oxygen.getLockingMask());
	CPPUNIT_ASSERT(MolpherAtom::NO_MUTATION & oxygen.getLockingMask());
	CPPUNIT_ASSERT(!(MolpherAtom::UNLOCKED & oxygen.getLockingMask()));
	CPPUNIT_ASSERT(!(MolpherAtom::KEEP_NEIGHBORS & oxygen.getLockingMask()));

	// copying
	MolpherAtom oxygen_copy = MolpherAtom(oxygen);
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getFormalCharge(), oxygen.getFormalCharge());
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getSymbol(), oxygen.getSymbol());
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getMass(), oxygen.getMass());
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getAtomicNum(), oxygen.getAtomicNum());
	CPPUNIT_ASSERT(oxygen_copy.isLocked());
	CPPUNIT_ASSERT(MolpherAtom::NO_ADDITION & oxygen_copy.getLockingMask());
	CPPUNIT_ASSERT(MolpherAtom::NO_MUTATION & oxygen_copy.getLockingMask());
	CPPUNIT_ASSERT(!(MolpherAtom::UNLOCKED & oxygen_copy.getLockingMask()));
	CPPUNIT_ASSERT(!(MolpherAtom::KEEP_NEIGHBORS & oxygen_copy.getLockingMask()));
	oxygen_copy = carbon;
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getFormalCharge(), carbon.getFormalCharge());
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getSymbol(), carbon.getSymbol());
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getMass(), carbon.getMass());
	CPPUNIT_ASSERT_EQUAL(oxygen_copy.getAtomicNum(), carbon.getAtomicNum());
	CPPUNIT_ASSERT(!carbon.isLocked());
	CPPUNIT_ASSERT(!(MolpherAtom::NO_ADDITION & oxygen_copy.getLockingMask()));
	CPPUNIT_ASSERT(!(MolpherAtom::NO_MUTATION & oxygen_copy.getLockingMask()));
	CPPUNIT_ASSERT(!(MolpherAtom::UNLOCKED & oxygen_copy.getLockingMask()));
	CPPUNIT_ASSERT(!(MolpherAtom::KEEP_NEIGHBORS & oxygen_copy.getLockingMask()));

}

void MinimalTest::testMolpherMol() {
    // test the simple constructor
    MolpherMol mol("CCO");
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(DBL_MAX, mol.getDistToTarget());

    // test the copy method
    auto mol_copy = mol.copy();
    mol_copy->setDistToTarget(1.0);
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol_copy->getSMILES());
    CPPUNIT_ASSERT_EQUAL(1.0, mol_copy->getDistToTarget());
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(DBL_MAX, mol.getDistToTarget());

    // molecules created by themselves should not belong to a tree
    CPPUNIT_ASSERT(!mol.getTree());
    CPPUNIT_ASSERT(!mol_copy->getTree());
    CPPUNIT_ASSERT(!mol.isBoundToTree());
    CPPUNIT_ASSERT(!mol_copy->isBoundToTree());

    mol.setSMILES("C1CCC1");
    CPPUNIT_ASSERT_EQUAL(std::string("C1CCC1"), mol.getSMILES());
    CPPUNIT_ASSERT_EQUAL(std::string("CCO"), mol_copy->getSMILES());

//    CPPUNIT_ASSERT_THROW(mol.setSMILES("ABCDE");, RDKit::SmilesParseException); // rdkit also gives a segmentation fault for some reason...

    // test the empty constructor (molecule with empty SMILES should not be valid and all values should be initialized to defaults)
    MolpherMol empty;
    CPPUNIT_ASSERT_EQUAL(std::string(""), empty.getSMILES());
    CPPUNIT_ASSERT(!empty.isValid());
    CPPUNIT_ASSERT_EQUAL(std::string(""), empty.getParentOper());
    CPPUNIT_ASSERT_EQUAL((unsigned) 0, empty.getItersWithoutDistImprovement());

    // test initialization from SDF
    RDKit::RWMol* mol_rdkit = new RDKit::RWMol(*(RDKit::SDMolSupplier(test_dir + "Structure2D_CID_4914.sdf").next()));
    RDKit::MolOps::Kekulize(*mol_rdkit);
    MolpherMol sdf_derived(test_dir + "Structure2D_CID_4914.sdf");
    CPPUNIT_ASSERT_EQUAL(std::string("CCN(CC)CCOC(=O)C1=CC=C(N)C=C1"), sdf_derived.getSMILES());

    std::ifstream t(test_dir + "Structure2D_CID_4914.sdf");
    std::string file_string(
            (std::istreambuf_iterator<char>(t)),
            std::istreambuf_iterator<char>());
    MolpherMol stream_derived(file_string);
    CPPUNIT_ASSERT_EQUAL(std::string("CCN(CC)CCOC(=O)C1=CC=C(N)C=C1"), stream_derived.getSMILES());

    delete mol_rdkit;

    // atom access and locking features
    RDKit::ROMol* mol_rdro = RDKit::SDMolSupplier(test_dir + "Structure2D_CID_4914.sdf").next();
    RDKit::RWMol* mol_rdrw = new RDKit::RWMol(*mol_rdro);
    MolpherMol mol_mphr_ro(mol_rdro);
    MolpherMol mol_mphr_rw(mol_rdrw);
    CPPUNIT_ASSERT(!mol_rdrw);

    int idx = 0;
    RDKit::RWMol* mol_rd_copy = mol_mphr_rw.asRDMol();
    for (auto atom : mol_mphr_rw.getAtoms()) {
        CPPUNIT_ASSERT_EQUAL(atom->getSymbol(), mol_rd_copy->getAtomWithIdx(idx)->getSymbol());
        CPPUNIT_ASSERT_EQUAL(atom->getSymbol(), mol_mphr_ro.getAtom(idx)->getSymbol());
        CPPUNIT_ASSERT_EQUAL(atom->getSymbol(), mol_rdro->getAtomWithIdx(idx)->getSymbol());

        CPPUNIT_ASSERT_EQUAL(atom->getLockingMask(), mol_mphr_ro.getAtom(idx)->getLockingMask());
        CPPUNIT_ASSERT(atom->getLockingMask() == MolpherAtom::UNLOCKED);

        idx++;
    }

    mol_mphr_rw.lockAtom(5, MolpherAtom::NO_MUTATION | MolpherAtom::KEEP_NEIGHBORS);
    mol_mphr_rw.getAtom(6)->setLockingMask(MolpherAtom::NO_ADDITION);
    idx = 0;
    for (auto atom : mol_mphr_rw.getAtoms()) {
        if (idx == 5) {
            CPPUNIT_ASSERT(atom->getLockingMask() == (MolpherAtom::NO_MUTATION | MolpherAtom::KEEP_NEIGHBORS | MolpherAtom::NO_REMOVAL));
        }
        if (idx == 6) {
            CPPUNIT_ASSERT(atom->getLockingMask() == MolpherAtom::NO_ADDITION);
        }
        idx++;
    }
    CPPUNIT_ASSERT(mol_mphr_rw.getAtom(6)->getLockingMask() == MolpherAtom::NO_ADDITION);
    CPPUNIT_ASSERT(mol_mphr_rw.getAtom(5)->getLockingMask() == (MolpherAtom::NO_MUTATION | MolpherAtom::KEEP_NEIGHBORS | MolpherAtom::NO_REMOVAL));

	// locks should be in the mol block
    std::string mol_block(mol_mphr_rw.asMolBlock(true));
    std::string expected_info_6(">  <MOLPHER_KEEP_NEIGHBORS>\n6\n");
	std::string expected_info_7(">  <MOLPHER_NO_ADDITION>\n7\n");
    CPPUNIT_ASSERT(mol_block.find(expected_info_6) != std::string::npos);
	CPPUNIT_ASSERT(mol_block.find(expected_info_7) != std::string::npos);

    MolpherMol ethanol_no_add(test_dir + "ethanol.sdf");
    std::set<int> locked_indices{1,2};
    std::vector<std::pair<int, std::shared_ptr<MolpherAtom>>> locked_atoms;
    idx = 0;
    for (auto atom : ethanol_no_add.getAtoms()) {
        if (atom->isLocked()) {
            locked_atoms.push_back(std::pair<int, std::shared_ptr<MolpherAtom>>(idx, atom));
        }
        idx++;
    }
    for (auto locked_atom : locked_atoms) {
        CPPUNIT_ASSERT_EQUAL(locked_atom.second->getSymbol(), std::string("C"));
        CPPUNIT_ASSERT(locked_indices.find(locked_atom.first) != locked_indices.end());
    }

    // test locking with smarts patterns and multiple indices
	MolpherMol malonic_acid("OC(=O)CC(=O)O");
    auto locked = malonic_acid.lockAtoms("OC(=O)", MolpherAtom::NO_ADDITION | MolpherAtom::NO_MUTATION);
	CPPUNIT_ASSERT(locked.size() == 6);
    for (auto i : locked) {
    	CPPUNIT_ASSERT_EQUAL(MolpherAtom::NO_ADDITION | MolpherAtom::NO_MUTATION, malonic_acid.getAtom(i)->getLockingMask());
    }
    locked = malonic_acid.lockAtoms("[$(*C(=O)O)]", MolpherAtom::NO_REMOVAL);
	CPPUNIT_ASSERT(locked.size() == 1);
	CPPUNIT_ASSERT_EQUAL((int) MolpherAtom::NO_REMOVAL, malonic_acid.getAtom(locked[0])->getLockingMask());
}

void MinimalTest::testAtomLibrary() {
    const AtomLibrary& default_lib = AtomLibrary::getDefaultLibrary();
	CPPUNIT_ASSERT_EQUAL(default_lib.getAtoms().size(), default_lib.getAtomProbabilities().size());
    std::vector<double> probabilities = default_lib.getAtomProbabilities();

    print("Default library: ");
    int idx = 0;
    for (auto atom : default_lib.getAtoms()) {
        RDKit::Atom* rd_atom = atom->asRDAtom();
        print(atom->getSymbol() + ", formal charge: " + std::to_string(atom->getFormalCharge())
              + ", probability: " + std::to_string(probabilities[idx]));
        idx += 1;
        CPPUNIT_ASSERT_EQUAL(atom->getSymbol(), rd_atom->getSymbol());
        CPPUNIT_ASSERT_EQUAL(atom->getMass(), rd_atom->getMass());
        CPPUNIT_ASSERT_EQUAL(atom->getAtomicNum(), (unsigned) rd_atom->getAtomicNum());
        CPPUNIT_ASSERT_EQUAL(atom->getFormalCharge(), rd_atom->getFormalCharge());
    }

    // changing the default library
	print("Modifying defaults... ");
    auto locked_nitrogen = std::make_shared<MolpherAtom>(MolpherAtom("N"));
    locked_nitrogen->setLockingMask(MolpherAtom::NO_ADDITION);
	AtomLibrary new_lib((std::vector<std::shared_ptr<MolpherAtom>>(
			{
					std::make_shared<MolpherAtom>(MolpherAtom("C"))
					, std::make_shared<MolpherAtom>(MolpherAtom("O"))
					, locked_nitrogen
					, std::make_shared<MolpherAtom>(MolpherAtom("F"))
					, std::make_shared<MolpherAtom>(MolpherAtom("S"))
			})), std::vector<double>(
			{
					73.12,
					11.741,
					11.318,
					1.379,
					1.295,
			})
	);
	AtomLibrary::setDefaultLibrary(new_lib);
	CPPUNIT_ASSERT_EQUAL(5, (int) default_lib.getAtomProbabilities().size());
	CPPUNIT_ASSERT_EQUAL(default_lib.getAtoms().size(), default_lib.getAtomProbabilities().size());

    print("New library: ");
	probabilities = default_lib.getAtomProbabilities();
    for (auto atom : default_lib.getAtoms()) {
        RDKit::Atom* rd_atom = atom->asRDAtom();
        if (default_lib.getAtomProbabilities().empty()) {
            print(atom->getSymbol() + ", formal charge: " + std::to_string(atom->getFormalCharge()));
        }else{
            print(atom->getSymbol() + ", formal charge: " + std::to_string(atom->getFormalCharge())
            + ", atom probabilities: " + std::to_string(probabilities[idx]));
            idx += 1;
        }
        CPPUNIT_ASSERT_EQUAL(atom->getSymbol(), rd_atom->getSymbol());
        CPPUNIT_ASSERT_EQUAL(atom->getMass(), rd_atom->getMass());
        CPPUNIT_ASSERT_EQUAL(atom->getAtomicNum(), (unsigned) rd_atom->getAtomicNum());
        CPPUNIT_ASSERT_EQUAL(atom->getFormalCharge(), rd_atom->getFormalCharge());
        if (atom->getSymbol() == "N") {
            CPPUNIT_ASSERT(MolpherAtom::NO_ADDITION & atom->getLockingMask());
            CPPUNIT_ASSERT(atom->isLocked());
        }
    }

    // getting a random atom
    print("Choosing random atoms until we generate the expected locked nitrogen: ");
    MolpherAtom rand_atom = default_lib.getRandomAtom();
    while (!rand_atom.isLocked()) {
        print(rand_atom.getSymbol());
        rand_atom = default_lib.getRandomAtom();
    }
    print("Locked nitrogen generated: " + rand_atom.getSymbol());
    CPPUNIT_ASSERT_EQUAL(rand_atom.getSymbol(), locked_nitrogen->getSymbol());
    CPPUNIT_ASSERT(MolpherAtom::NO_ADDITION & locked_nitrogen->getLockingMask());
    CPPUNIT_ASSERT(rand_atom.isLocked());
}

void MinimalTest::testAddAtomOperator() {
    std::shared_ptr<MolpherMol> cymene_no_add(new MolpherMol(test_dir + "cymene.sdf"));
	print("original\n" + cymene_no_add->getSMILES() + "\n");
	print_lock_info(cymene_no_add);

    AddAtom op_add;
    op_add.setOriginal(cymene_no_add);
    const std::vector<unsigned int> open_atom_indices = op_add.getOpenIndices();
	std::vector<std::shared_ptr<MolpherAtom>> open_atoms = op_add.getOpenAtoms();
	int counter = 0;
	for (auto idx : open_atom_indices) {
		MolpherAtom* atom_orig = cymene_no_add->getAtom(idx).get();
		MolpherAtom* atom_ret = open_atoms[counter].get();
		CPPUNIT_ASSERT_EQUAL(atom_orig, atom_ret);
		counter++;
	}

    // first gen
	auto first_morph = op_add.morph();
	int idx = 0;
	for (auto atom : cymene_no_add->getAtoms()) {
		CPPUNIT_ASSERT_EQUAL(first_morph->getAtom(idx)->getLockingMask(), atom->getLockingMask());
		idx++;
	}
	counter = 200;
    while (counter != 0) {
        CPPUNIT_ASSERT(match_substr(op_add.morph()->getSMILES(), "c1ccccc1"));
        counter--;
    }

    // second gen
	op_add.setOriginal(first_morph);
    counter = 200;
    while (counter != 0) {
        CPPUNIT_ASSERT(match_substr(op_add.morph()->getSMILES(), "c1ccccc1"));
        counter--;
    }

	// custom library
	AtomLibrary custom_lib((std::vector<std::shared_ptr<MolpherAtom>>(
			{
					std::make_shared<MolpherAtom>(MolpherAtom("O"))
					, std::make_shared<MolpherAtom>(MolpherAtom("S"))
			}))
	);
	op_add = AddAtom(custom_lib);
	op_add.setOriginal(cymene_no_add);
	counter = 200;
	while (counter != 0) {
		auto morph = op_add.morph();
		print(morph->getSMILES());
		CPPUNIT_ASSERT(match_substr(morph->getSMILES(), "c1ccccc1"));
		CPPUNIT_ASSERT(morph->getSMILES().find("S") != std::string::npos || morph->getSMILES().find("O") != std::string::npos);
		counter--;
	}
}

void MinimalTest::testRemoveAtomOperator() {
	std::shared_ptr<MolpherMol> dimethylaniline(new MolpherMol(test_dir + "dimethylaniline.sdf"));
	print("original\n" + dimethylaniline->getSMILES() + "\n");
	print_lock_info(dimethylaniline);

	CPPUNIT_ASSERT_EQUAL(dimethylaniline->getAtom(3)->getLockingMask(), MolpherAtom::KEEP_NEIGHBORS | MolpherAtom::NO_REMOVAL);
	CPPUNIT_ASSERT_EQUAL(dimethylaniline->getAtom(8)->getLockingMask(), (int) MolpherAtom::UNLOCKED);

	RemoveAtom op_remove;
	op_remove.setOriginal(dimethylaniline);
	print("original\n" + dimethylaniline->getSMILES() + "\n");
	for (int i = 0; i != 10; i++) {
		CPPUNIT_ASSERT_EQUAL(std::string("CC1=C(N)C=CC=C1"), op_remove.morph()->getSMILES());
	}

	auto morph = op_remove.morph();
	op_remove.setOriginal(morph);
	CPPUNIT_ASSERT(op_remove.getMarkedIndices().empty());
	CPPUNIT_ASSERT(op_remove.getMarkedAtoms().empty());
	for (int i = 0; i != 10; i++) {
		CPPUNIT_ASSERT(!op_remove.morph());
	}

	std::shared_ptr<MolpherMol> dimethylphenylformamide(new MolpherMol(test_dir + "dimethylphenylformamide.sdf"));
	op_remove.setOriginal(dimethylphenylformamide);
	CPPUNIT_ASSERT_EQUAL(dimethylphenylformamide.get(), op_remove.getOriginal().get());
	RDKit::ROMol *mol1 = nullptr;
	RDKit::RWMol *patt = nullptr;
	for (int i = 0; i != 10; i++) {
		morph = op_remove.morph();
		if (!morph) {
			op_remove.setOriginal(dimethylphenylformamide);
		} else {
			print(morph->getSMILES());
			std::string smiles = morph->getSMILES();
			mol1 = RDKit::SmilesToMol(smiles);
			patt = RDKit::SmartsToMol("CcccC");
			RDKit::MatchVectType res;
			CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt , res ));
			op_remove.setOriginal(morph);
		}
	}

	delete mol1;
	delete patt;
}

void MinimalTest::testAddBondOperator() {
	std::shared_ptr<MolpherMol> propanol(new MolpherMol(test_dir + "propanol.sdf"));
	print("original\n" + propanol->getSMILES() + "\n");
	print_lock_info(propanol);

	unsigned int locked_idx = 2;
	CPPUNIT_ASSERT_EQUAL(propanol->getAtom(locked_idx)->getLockingMask(), (int) (MolpherAtom::KEEP_NEIGHBORS_AND_BONDS | MolpherAtom::KEEP_NEIGHBORS | MolpherAtom::KEEP_BONDS | MolpherAtom::NO_REMOVAL));

	AddBond op_add_bond;
	op_add_bond.setOriginal(propanol);
	print("original\n" + propanol->getSMILES() + "\n");
	for (auto bond : op_add_bond.getOpenBonds()) {
		CPPUNIT_ASSERT(bond.first != locked_idx && bond.second != locked_idx);
	}

	RDKit::ROMol* propanol_rd = propanol->asRDMol();
	for (int i = 0; i != 10; i++) {
		auto morph = op_add_bond.morph();
		print(morph->getSMILES());
		RDKit::ROMol* morph_rd = morph->asRDMol();

		for (int atom_idx = 0; atom_idx != propanol->getAtomCount(); atom_idx++) {
			CPPUNIT_ASSERT_EQUAL(propanol->getAtom(atom_idx)->getLockingMask(), morph->getAtom(atom_idx)->getLockingMask());
			if (atom_idx == locked_idx) {
				CPPUNIT_ASSERT_EQUAL(propanol_rd->getAtomWithIdx(atom_idx)->getNumImplicitHs(), morph_rd->getAtomWithIdx(atom_idx)->getNumImplicitHs());
			}
		}
		delete morph_rd;
	}
	delete propanol_rd;
}

void MinimalTest::testRemoveBondOperator() {
	std::shared_ptr<MolpherMol> test_mol(new MolpherMol(test_dir + "remove_bond_test_mol.sdf"));
	print("original\n" + test_mol->getSMILES() + "\n");
	print_lock_info(test_mol);

	RemoveBond op_remove_bond;
	op_remove_bond.setOriginal(test_mol);

	RDKit::ROMol* test_mol_rd = test_mol->asRDMol();
	RDKit::ROMol *mol1 = nullptr;
	RDKit::RWMol *patt = nullptr;
	for (int i = 0; i != 20; i++) {
		auto morph = op_remove_bond.morph();
		if (morph == nullptr) break;
		print(morph->getSMILES());
		RDKit::ROMol* morph_rd = morph->asRDMol();

		std::string smiles = morph->getSMILES();
		mol1 = RDKit::SmilesToMol(smiles);
		patt = RDKit::SmartsToMol("c1ccccc1");
		RDKit::MatchVectType res;
		CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt , res ));
//		op_remove_bond.setOriginal(morph);
		delete morph_rd;
	}

	delete mol1;
	delete patt;
	delete test_mol_rd;
}

void MinimalTest::testMutateAtomOperator() {
	std::shared_ptr<MolpherMol> alanine(new MolpherMol(test_dir + "alanine.sdf"));
	print("original\n" + alanine->getSMILES() + "\n");
	print_lock_info(alanine);

	MutateAtom op_mutate;
	op_mutate.setOriginal(alanine);
	// TODO: expand this test

	RDKit::ROMol* test_mol_rd = alanine->asRDMol();
	RDKit::ROMol *mol1 = nullptr;
	RDKit::RWMol *patt = nullptr;
	for (int i = 0; i != 50; i++) {
		auto morph = op_mutate.morph();
		if (morph == nullptr) continue;
		print(morph->getSMILES());
		RDKit::ROMol* morph_rd = morph->asRDMol();

		mol1 = RDKit::SmilesToMol(morph->getSMILES());
		patt = RDKit::SmartsToMol("C(=O)");
		RDKit::MatchVectType res;
		CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt , res ));
//		op_remove_bond.setOriginal(morph);
		delete morph_rd;
	}

	delete mol1;
	delete patt;
	delete test_mol_rd;
}

void MinimalTest::testInterlayAtomOperator() {
	std::shared_ptr<MolpherMol> isopropylphenol(new MolpherMol(test_dir + "isopropylphenol.sdf"));
	print("original\n" + isopropylphenol->getSMILES() + "\n");
	print_lock_info(isopropylphenol);

	InterlayAtom op_interlay;
	op_interlay.setOriginal(isopropylphenol);
	// TODO: expand this test

	RDKit::ROMol* test_mol_rd = isopropylphenol->asRDMol();
	RDKit::ROMol *mol1 = nullptr;
	RDKit::RWMol *patt = nullptr;
	for (int i = 0; i != 50; i++) {
		auto morph = op_interlay.morph();
		if (morph == nullptr) continue;
		print(morph->getSMILES());
		RDKit::ROMol* morph_rd = morph->asRDMol();

		mol1 = RDKit::SmilesToMol(morph->getSMILES());
		patt = RDKit::SmartsToMol("C(C)C");
		RDKit::MatchVectType res;
		CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt , res ));
//		op_interlay.setOriginal(morph);
		delete morph_rd;
	}

	delete mol1;
	delete patt;
	delete test_mol_rd;
}

void MinimalTest::testContractBondOperator() {
	std::shared_ptr<MolpherMol> test_mol(new MolpherMol(test_dir + "contract_bond_test_mol.sdf"));
	print("original\n" + test_mol->getSMILES() + "\n");
	print_lock_info(test_mol);

	ContractBond op_contract;
	op_contract.setOriginal(test_mol);
	// TODO: expand this test (add a molecule with a double bond)

	RDKit::ROMol* test_mol_rd = test_mol->asRDMol();
	for (int i = 0; i != 50; i++) {
		auto morph = op_contract.morph();
		if (morph == nullptr) continue;
		print(morph->getSMILES());

		CPPUNIT_ASSERT(morph->getSMILES() == "CCCC1CC(N)CC1C(=O)O" || morph->getSMILES() == "CCCC1CC(C(=O)O)CC1N");
	}

	print("lock test");
	auto morph = op_contract.morph();
	for (auto atm : morph->getAtoms()) {
		std::cout << atm->getSymbol() << ": " << atm->getLockingMask() << std::endl;
	}
	print(morph->getSMILES());
	op_contract.setOriginal(morph);
	morph = op_contract.morph();
	print(morph->getSMILES());
	CPPUNIT_ASSERT ("CCCC1C(N)CC1C(=O)O" == morph->getSMILES());

	op_contract.setOriginal(morph);
	morph = op_contract.morph();
	CPPUNIT_ASSERT(morph == nullptr);

	delete test_mol_rd;
}

void MinimalTest::testRerouteBondOperator() {
	std::shared_ptr<MolpherMol> test_mol(new MolpherMol(test_dir + "reroute_test.sdf"));
	print("original\n" + test_mol->getSMILES() + "\n");
	print_lock_info(test_mol);

	RerouteBond op_reroute;
	op_reroute.setOriginal(test_mol);

	RDKit::ROMol *mol1 = nullptr;
	RDKit::RWMol *patt = nullptr;
	for (int i = 0; i != 50; i++) {
		auto morph = op_reroute.morph();
		if (morph == nullptr) continue;
		print(morph->getSMILES());

		RDKit::ROMol* morph_rd = morph->asRDMol();

		mol1 = RDKit::SmilesToMol(morph->getSMILES());
		patt = RDKit::SmartsToMol("C(F)C(O)C(C)Cl");
		RDKit::MatchVectType res;
		CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt , res ));
		op_reroute.setOriginal(morph);
		delete morph_rd;
	}

	delete mol1;
	delete patt;
}

void MinimalTest::testMolpher() {
	std::vector<std::shared_ptr<MorphingOperator>> opers = {
			std::make_shared<AddAtom>()
			, std::make_shared<RemoveAtom>()
			, std::make_shared<AddBond>()
			, std::make_shared<RemoveBond>()
			, std::make_shared<MutateAtom>()
			, std::make_shared<InterlayAtom>()
			, std::make_shared<ContractBond>()
			, std::make_shared<RerouteBond>()
	};
	std::shared_ptr<MolpherMol> dimethylphenylformamide(new MolpherMol(test_dir + "dimethylphenylformamide_complex_locking.sdf"));
	print("original\n" + dimethylphenylformamide->getSMILES() + "\n");
	print_lock_info(dimethylphenylformamide);

	Molpher molpher(dimethylphenylformamide, opers, 0, 200);
	molpher();
	RDKit::ROMol *mol1 = nullptr;
	RDKit::RWMol *patt2 = nullptr;
	for (auto morph : molpher.getMorphs()) {
		print(morph->getSMILES());
		mol1 = RDKit::SmilesToMol(morph->getSMILES());
		patt2 = RDKit::SmartsToMol("c1ccccc1");
		RDKit::MatchVectType res;
		CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt2 , res ));
	}
	CPPUNIT_ASSERT(molpher.getMorphs().empty());
	CPPUNIT_ASSERT_EQUAL(dimethylphenylformamide.get(), molpher.getOriginal().get());

	std::shared_ptr<MolpherMol> cymene(new MolpherMol(test_dir + "cymene.sdf"));
	molpher.reset(cymene);
	CPPUNIT_ASSERT_EQUAL(cymene.get(), molpher.getOriginal().get());
	print("original\n" + cymene->getSMILES() + "\n");
	print_lock_info(cymene);

	molpher();
	for (auto morph : molpher.getMorphs()) {
		print(morph->getSMILES());
		mol1 = RDKit::SmilesToMol(morph->getSMILES());
		patt2 = RDKit::SmartsToMol("c1ccccc1");
		RDKit::MatchVectType res;
		CPPUNIT_ASSERT(RDKit::SubstructMatch( *mol1 , *patt2 , res ));
	}

	// test use of custom operator
	opers.push_back(std::make_shared<IdentityOperator>());
	Molpher molpher_custom_oper(cymene, opers, 0, 200);
	molpher_custom_oper();
	bool found_orig(false);
	for (auto morph : molpher_custom_oper.getMorphs()) {
		if (morph->getSMILES() == cymene->getSMILES()) {
			found_orig = true;
		}
	}
	CPPUNIT_ASSERT(found_orig);

	// test collector interface
	std::vector<std::shared_ptr<MorphCollector> > collectors;
	collectors.push_back(std::static_pointer_cast<MorphCollector>(std::make_shared<PrintMorphingInfo>()));
	Molpher molpher_custom_oper_collectors(cymene, opers, 0, 200, collectors);
	molpher_custom_oper_collectors();
	for (auto morph : molpher_custom_oper_collectors.getMorphs()) {
		CPPUNIT_ASSERT_EQUAL(cymene->getSMILES(), morph->getParentSMILES());
		CPPUNIT_ASSERT_EQUAL(morph->getSAScore(), morph->getDistToTarget());
	}

	delete mol1;
	delete patt2;
}

void MinimalTest::testExplorationData() {
    // test empty initialzation, default data should not be valid
    ExplorationData data;
    CPPUNIT_ASSERT(!data.isValid());
    
    // get the source molecule, should be empty
    CPPUNIT_ASSERT(!data.getSource());
	CPPUNIT_ASSERT(!data.getTarget());
    
    // modify the source and target molecule and save them back to tree data
	auto source(std::make_shared<MolpherMol>("NC(=O)C"));
	auto target(std::make_shared<MolpherMol>("CC(=O)C"));
    data.setSource(*source);
    data.setTarget(*target);
    
    // data should now be valid
    CPPUNIT_ASSERT(data.isValid());
    
    // the tree should now contain the source as root
    auto map = data.getTreeMap();
    CPPUNIT_ASSERT(map->size() == 1);
    CPPUNIT_ASSERT_EQUAL(map->find(source->getSMILES())->second->getSMILES(), source->getSMILES());
    
    // we shouldn't be able to set source molecule for the second time
    CPPUNIT_ASSERT_THROW(data.setSource(*source);, std::runtime_error);
    data.setTarget(*target); // we can change the target though
    
    // check the selectors
    CPPUNIT_ASSERT_EQUAL(FP_MORGAN, static_cast<FingerprintSelector>(data.getFingerprint()));
    CPPUNIT_ASSERT_EQUAL(SC_TANIMOTO, static_cast<SimCoeffSelector>(data.getSimilarityCoefficient()));
    auto operators = data.getChemicalOperators();
    CPPUNIT_ASSERT(!operators.empty());
    for (auto& oper : operators) {
        std::cout << ChemOperLongDesc(oper) << std::endl;
    }
    
    // test the serialization
    auto data_from_template = ExplorationData::load(test_dir + "test-template.xml");
    CPPUNIT_ASSERT(data_from_template->isValid());
    
    data_from_template->setTarget(*target);
    data_from_template->save(test_dir + "template_snapshot.xml");
    CPPUNIT_ASSERT_THROW(ExplorationData::load("./does/not/exist"), std::runtime_error);
    auto data_from_snapshot = ExplorationData::load(test_dir + "template_snapshot.xml");
    CPPUNIT_ASSERT_EQUAL(target->getSMILES(), data_from_snapshot->getTarget()->getSMILES());
}


void MinimalTest::testTree() {
    // tree creation
    auto tree = ExplorationTree::create("CCO", "C1=COC=C1");
    CPPUNIT_ASSERT_THROW(ExplorationTree::create("CCO", "");, std::runtime_error);
    
    // retrieve the source and see if it belongs to the correct tree
	CPPUNIT_ASSERT_THROW(tree->fetchMol("CC");, std::runtime_error);
    auto source = tree->fetchMol("CCO");
    CPPUNIT_ASSERT(source->isBoundToTree());
    auto tree_of_source = source->getTree();
    CPPUNIT_ASSERT(tree_of_source);
    CPPUNIT_ASSERT_EQUAL(tree, tree_of_source);
    CPPUNIT_ASSERT(tree_of_source->hasMol(source));
    CPPUNIT_ASSERT(tree_of_source->hasMol(source->getSMILES()));
    
    // we should not be able to set ownership of an already owned molecule
    CPPUNIT_ASSERT_THROW(source->setOwner(tree);, std::runtime_error);
    
    // fetch leaves shouldn't increase the distance imporvement counter for source
    auto leaves = tree->fetchLeaves(true);
    CPPUNIT_ASSERT(leaves.size() == 1);
    CPPUNIT_ASSERT_EQUAL(source, leaves[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned) 0, source->getItersWithoutDistImprovement());
    
    // generate some candidate morphs
    CPPUNIT_ASSERT(tree->getCandidateMorphs().empty());
    CPPUNIT_ASSERT(tree->getCandidateMorphsMask().empty());
    tree->generateMorphs();
    CPPUNIT_ASSERT(!tree->getCandidateMorphs().empty());
    CPPUNIT_ASSERT(!tree->getCandidateMorphsMask().empty());
    CPPUNIT_ASSERT_EQUAL(tree->getCandidateMorphs().size(), tree->getCandidateMorphsMask().size());
    for (auto candidate : tree->getCandidateMorphs()) {
        CPPUNIT_ASSERT(!candidate->isBoundToTree());
        CPPUNIT_ASSERT_EQUAL(false, (bool) candidate->getTree());
    }
    printCandidates(tree);
    
    // sort the morphs
    tree->sortMorphs();
    printCandidates(tree);
    double previous = 0;
    for (auto candidate : tree->getCandidateMorphs()) {
        CPPUNIT_ASSERT(candidate->getDistToTarget() >= previous);
        previous = candidate->getDistToTarget();
    }
    
    // filter morphs
    tree->filterMorphs(true);
    printCandidates(tree);
    
    // extend the tree
    CPPUNIT_ASSERT_EQUAL((unsigned) 0, tree->getGenerationCount());
    
    tree->extend();
    leaves = tree->fetchLeaves();
    auto source_descendents = source->getDescendants();
    for (auto leaf : leaves) {
        CPPUNIT_ASSERT(leaf->isBoundToTree());
        CPPUNIT_ASSERT_EQUAL(tree, leaf->getTree());
        
        auto parent_smiles = leaf->getParentSMILES();
        CPPUNIT_ASSERT_EQUAL(source->getSMILES(), parent_smiles);
        CPPUNIT_ASSERT_EQUAL(source, tree->fetchMol(parent_smiles));
        
        auto it = source_descendents.find(leaf->getSMILES());
        CPPUNIT_ASSERT(it != source_descendents.end());
        CPPUNIT_ASSERT_EQUAL(leaf, tree->fetchMol(*it));
    }
    
    CPPUNIT_ASSERT(tree->getCandidateMorphs().empty());
    CPPUNIT_ASSERT(tree->getCandidateMorphsMask().empty());
    CPPUNIT_ASSERT_EQUAL((unsigned) 1, tree->getGenerationCount());
    
    // prune the tree
    tree->prune();
    auto leaf_zero = leaves[0];
    tree->deleteSubtree(leaf_zero->getSMILES());
    CPPUNIT_ASSERT(!tree->hasMol(leaf_zero->getSMILES()));
    CPPUNIT_ASSERT(!leaf_zero->getTree());
    leaf_zero->removeFromTree(); // shouldn't do anything
    
    // remove all newly generated morphs
    CPPUNIT_ASSERT_THROW(tree->deleteSubtree(source->getSMILES(), false);, std::runtime_error); // should not be possible
    tree->deleteSubtree(source->getSMILES(), true); // remove only descendents
    
    source_descendents = source->getDescendants();
    auto source_hist_descendents = source->getHistoricDescendants();
    for (auto leaf : leaves) {
        CPPUNIT_ASSERT(!tree->hasMol(leaf->getSMILES()));
        CPPUNIT_ASSERT(!leaf->getTree());
        
        auto it = source_descendents.find(leaf->getSMILES());
        CPPUNIT_ASSERT(it == source_descendents.end());
        CPPUNIT_ASSERT(it != source_hist_descendents.end());
    }
    
    // make a few more generations
    for (unsigned iter_idx = 0; iter_idx != 5; iter_idx++) {
        tree->generateMorphs();
        tree->sortMorphs();
        tree->filterMorphs();
        printCandidates(tree);
        tree->extend();
        std::cout << "Path found: " + NumberToStr(tree->isPathFound()) << std::endl;
        tree->prune();
    }
    
    // test the tree traversal
    PrintMols printing_callback;
    tree->traverse(printing_callback);
    
    // serialize the tree into file
    tree->save(test_dir + "testTree_snapshot.xml");
    
    // create new tree from the file
    auto tree_from_file = ExplorationTree::create(test_dir + "testTree_snapshot.xml");
    
    // find the path
    while (true) {
		tree_from_file->generateMorphs();
        tree_from_file->sortMorphs();
        tree_from_file->filterMorphs();
        printCandidates(tree_from_file);
        tree_from_file->extend();
        std::cout << "Path found: " + NumberToStr(tree_from_file->isPathFound()) << std::endl;
		if (tree_from_file->isPathFound()) break;
		tree_from_file->prune();
    }
}

void MinimalTest::testTreeOperatorsAndLocks() {
	auto source = std::make_shared<MolpherMol>(test_dir + "gonane.sdf");
	print(source->getSMILES());
	auto tree = ExplorationTree::create(source);

	// check if locks maintained in the tree
	auto source_from_tree(tree->fetchMol(source->getSMILES()));
	print_lock_info(source_from_tree);
	auto atoms = source_from_tree->getAtoms();
	for (int idx = 0; idx != atoms.size(); idx++) {
		CPPUNIT_ASSERT_EQUAL(source->getAtom(idx)->getLockingMask(), atoms[idx]->getLockingMask());
	}

	// test adjusting tree operators
	auto orig_opers = tree->getMorphingOperators();
	std::vector<std::shared_ptr<MorphingOperator>> init_opers = {
			std::make_shared<AddAtom>()
	};
	tree->setMorphingOperators(init_opers);

	// run the initial phase (adding atoms only)
	int filters = FilterMorphsOper::MorphFilters::WEIGHT
				  | FilterMorphsOper::MorphFilters::SYNTHESIS
				  | FilterMorphsOper::MorphFilters::DUPLICATES;
	for (int iter_idx = 0; iter_idx != 2; iter_idx++) {
		tree->generateMorphs();
		tree->filterMorphs((FilterMorphsOper::MorphFilters) filters, false);
		tree->extend();
		CPPUNIT_ASSERT(match_substr(tree->getCandidateMorphs(), source->getSMILES()));
	}

	// add back original operators and repeat, use a collector too this time
	tree->setMorphingOperators(orig_opers);
	std::vector<std::shared_ptr<MorphCollector> > collectors;
	collectors.push_back(std::static_pointer_cast<MorphCollector>(std::make_shared<PrintMorphingInfo>()));
	for (int iter_idx = 0; iter_idx != 2; iter_idx++) {
		tree->generateMorphs(collectors);
		for (auto morph : tree->getCandidateMorphs()) {
			CPPUNIT_ASSERT_EQUAL(morph->getSAScore(), morph->getDistToTarget());
		}
		tree->filterMorphs((FilterMorphsOper::MorphFilters) filters, false);
		tree->extend();
		CPPUNIT_ASSERT(match_substr(tree->getCandidateMorphs(), source->getSMILES()));
	}
}

