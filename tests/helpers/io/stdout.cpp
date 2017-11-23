//
// Created by sichom on 6/6/17.
//

#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "stdout.hpp"

void print_mol_info(RDKit::ROMol *mol) {
    std::string smiles = RDKit::MolToSmiles(*mol);
    std::cout << smiles << std::endl;
    uint atom_count = mol->getNumAtoms(true);
    RDKit::Atom* atom = nullptr;
    for (uint atom_idx = 0; atom_idx != atom_count; atom_idx++) {
        atom = mol->getAtomWithIdx(atom_idx);
        std::cout << atom_idx << " - " + atom->getSymbol() << ": degree: " << atom->getDegree() << ", implicit Hs:" << atom->getNumImplicitHs() << std::endl;
    }
}

void print_morphs(const MolpherMol& parent, const std::vector<std::shared_ptr<MolpherMol>>& mols) {
    std::cout << "Morphs generated with fixed atoms (from: " << parent.getSMILES() << "): " << mols.size() << std::endl;
    for (auto morph : mols) {
        std::cout << morph->getSMILES() << std::endl;
    }
}

void print_locks(const MolpherMol& mol) {
    for (auto atom : mol.getAtoms()) {
        print(atom->getSymbol());
        for (auto lock : MolpherAtom::lockingMaskToString(atom->getLockingMask())) {
            print(lock);
        }
    }
}

void print(const std::string &text) {
    std::cout << text << std::endl;
}


