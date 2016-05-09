/*
 Copyright (c) 2013 Milan Vorsilak

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

//#include "main.hpp"
#include "SAScore.h"
#include "global_types.h"
#include "core/chem/fingerprintStrategy/FingerprintStrategy.h"
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>

SAScore * SAScore::instance = NULL;

/*
 Static method for loading data into instance.
 */
void SAScore::loadData() {
    SAScore::loadData("SAScore.dat");
}

void SAScore::loadData(const std::string& path) {
    std::cout << "loading SAScore.dat ... ";
    
    SAScore* inst = getInstance();
    std::ifstream myfile;
    // This file must be in the same directory as exe file of server.
    myfile.open(path.c_str());
    if (!myfile.good()) {
        // problem when opening the file
        throw std::runtime_error("Can't load sascore file SAScore.dat");
        return;
    }
    
    char array[31];
    myfile.getline(array, 31);

    unsigned int id;
    unsigned int count;
    double score;
    int i =0;
    while (!myfile.eof()) {
        myfile >> id >> count >> score;
        inst->data[id] = score;
        i++;
    }

    myfile.close();
    std::cout << "done" << std::endl;
    return;
}


/*
 This class is implemented as singleton.
 */
SAScore* SAScore::getInstance()
{
    if (instance == NULL) {
        instance = new SAScore();
    }    
    return instance;
}

SAScore::SAScore()
{
}

/*
 * Method returning sascore of given molecule like in http://www.jcheminf.com/content/1/1/8
 * Not optimized.
 */
double SAScore::getScore(RDKit::ROMol& mol)
{
    RDKit::SparseIntVect< boost::uint32_t >::StorageType::const_iterator iter;
    RDKit::SparseIntVect< boost::uint32_t > * fp = RDKit::MorganFingerprints::getFingerprint(mol, 2, 0, 0, false, true, false, 0);
    iter=fp->getNonzeroElements().begin();
    // compute fragment score
    double fragmentScore = 0;
    int sumOfFragments = 0;
    while (iter != fp->getNonzeroElements().end()) {
        //std::cout << iter->first << ":" << iter->second << "\n";
        if (data.find(iter->first) != data.end()) {
            //std::cout << iter->first << ": " << data[iter->first] << std::endl;
            fragmentScore+=data[iter->first]*iter->second;
            sumOfFragments+=iter->second;
        } else {
            fragmentScore-=2.95952;
        }
        iter++;
    }
    fragmentScore/=sumOfFragments;
    // compute complexity penalty
    double sizePenalty = pow(mol.getNumAtoms(true), 1.005) - mol.getNumAtoms(true);
    RDKit::VECT_INT_VECT rings = mol.getRingInfo()->atomRings();

    int macroCycleCount = 1;
    int bridgeAtomCount = 1;
    int spiroAtomCount = 1;
    int i=0;
    for (RDKit::VECT_INT_VECT_CI ringIter = rings.begin();ringIter != rings.end();ringIter++) {
        if (ringIter->size()>8) {
            macroCycleCount++;
        }
        //std::cout << "Cycle" << std::endl;
        //std::vector<int> i = rings[0].begin();
        //i = ringIter->begin();
        std::sort(rings[i].begin(),rings[i].end());
        for (RDKit::INT_VECT_CI it=ringIter->begin();it!=ringIter->end();it++) {
            //std::cout << *it << std::endl;
        }
        int x=1;
        for (RDKit::VECT_INT_VECT_CI ringIter2 = ringIter+1; ringIter2 != rings.end(); ringIter2++) {
            std::vector<int> v;
            RDKit::Intersect(rings[i],rings[i+x],v);
            if (v.size() == 1) {
                spiroAtomCount++;
            } else if (v.size() > 2) {
                bridgeAtomCount++;
            }
            x++;
        }
        i++;
    }
    // macrocycles are prohibited!
    //double macroCyclePenalty = log10(macroCycleCount);
    double ringComplexityPenalty = log10(bridgeAtomCount) + log10(spiroAtomCount);

    // release data
    delete fp;
    fp = 0;    
    
    if (macroCycleCount>1) {
        // some bigger value than max allowed score (6)
        
        // TODO: Use some well defined value here instead
        return 10;
    } else {
        return (sizePenalty+ringComplexityPenalty-fragmentScore) * 2 + 7.5;
    }
}

SAScore::~SAScore()
{
    data.clear();
}

SAScore* SAScore::destroyInstance()
{    
    if (instance != NULL) 
    {
        delete instance;
    }
    instance = NULL;
}

