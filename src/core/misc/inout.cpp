/*
 Copyright (c) 2012 Petr Koupy
 Copyright (c) 2012 Vladimir Fiklik
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

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <boost/thread/mutex.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/BadFileException.h>

#include "selectors/chemoper_selectors.h"
#include "core/misc/inout.h"

boost::mutex gIoMutex;

void SynchCout(const std::string &s) {
    boost::mutex::scoped_lock lock(gIoMutex);
    std::cout << s << std::endl;
}

void SynchCerr(const std::string &s, const std::string prefix) {
    boost::mutex::scoped_lock lock(gIoMutex);
    std::cerr << prefix << s << std::endl;
}

void Cout(const std::string &s) {
    std::cout << s << std::endl;
}

void Cerr(const std::string &s, const std::string prefix) {
    std::cerr << prefix << s << std::endl;
}

/**
 * Read molecules from txt file where each line contains just one smile
 * without white spaces.
 * @param file
 * @param mols
 */
void ReadRWMolFromTxtFile(const std::string &file,
        std::vector<RDKit::RWMol *> &mols) {
    std::ifstream stream;
    stream.open(file.c_str());

    if (stream.is_open()) {
		while (stream.good())  {
            std::string line;
			std::getline(stream, line);
            if (line.empty()) {
                // skip empty line
                continue;
            }
			// parse smiles ..
            try {
                RDKit::RWMol* mol = RDKit::SmilesToMol(line);
                // add to result
                mols.push_back(mol);
            } catch (std::exception e) {
                // failed to load the moll
                SynchCout("Failed to load mol from smile: " + line);
            }
		}
		stream.close();
	} else {
        // failed to open file
    }
}

void WriteRWMolsToSDF(const std::string &file,
        std::vector<RDKit::RWMol *> &mols) {
    try {
        RDKit::SDWriter writer(file);

        std::vector<RDKit::RWMol *>::iterator it;
        for (it = mols.begin(); it != mols.end(); ++it) {
            RDKit::RWMol *mol = *it;
            writer.write(*mol);
        }

        writer.close();
    } catch (RDKit::BadFileException &exc) {
        SynchCout(std::string("Cannot write to file: ").append(file));
    }
}

void ReadRWMolsFromSDF(const std::string &file,
        std::vector<RDKit::RWMol *> &mols) {
    try {
        RDKit::SDMolSupplier supplier(file);

        mols.reserve(supplier.length());
        for (size_t i = 0; i < supplier.length(); ++i) {
            RDKit::ROMol *mol = supplier.next();
            if (mol) {
                RDKit::RWMol *rwMol = new RDKit::RWMol(*mol);
                if (rwMol) {
                    mols.push_back(rwMol);
                }
                delete mol;
            }
        }
    } catch (RDKit::BadFileException &exc) {
        SynchCout(std::string("Cannot read from file: ").append(file));
    }
}
