//
// Created by sichom on 6/9/17.
//

#include <mol_helpers.hpp>
#include "io/stdout.hpp"
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

bool match_substr(std::string smiles_mol, std::string smiles_query) {
	RDKit::RWMol* mol(RDKit::SmilesToMol(smiles_mol));
	RDKit::RWMol* query(RDKit::SmilesToMol(smiles_query));

	RDKit::MatchVectType match;
	bool matches = RDKit::SubstructMatch(*mol, *query, match);

	return matches;
}

void print_lock_info(std::shared_ptr<MolpherMol> mol) {
	int counter = 0;
	for (auto atm : mol->getAtoms()) {
		std::vector<std::string> xx = MolpherAtom::lockingMaskToString(atm->getLockingMask());
		print(atm->getSymbol() + "." + std::to_string(++counter));
		for (auto x : xx) {
			print(x);
		}
	}
}
