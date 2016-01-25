/*
 * File:   VectorFpFngpr.cpp
 * Author: Petr Škoda
 *
 */

#include "VectorFpFngpr.hpp"
#include "global_types.h"

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

#include <boost/thread/mutex.hpp>
boost::mutex gFpIoMutex;

typedef std::vector<std::pair<boost::uint32_t,boost::uint32_t> > VectorInfo;

// Cache for fragments.
std::map<std::string, boost::shared_ptr<ExplicitBitVect> > cache;

VectorFp::VectorFp(unsigned int radius, unsigned int nChunks) : mRadius(radius), mNChunks(nChunks)
{
    // on-op
}

/**
 * Store given value in binary representation to given fingerprint.
 * @param value
 * @param size
 * @param fingerprint
 * @param offset
 */
void storeValueBinaryEncoding(unsigned int value, int size, ExplicitBitVect& fingerprint,
        int offset) {
    for (int i = 0; i < size; ++i) {
        if ( (value & (1 >> i)) != 0) {
            fingerprint.setBit(offset + i);
        }
    }
}

Fingerprint *VectorFp::GetFingerprint(RDKit::ROMol *mol) {    
    boost::mutex::scoped_lock lock(gFpIoMutex);

    // - - - - - VECTORFP SETTINGS - - - - -
    const int chunkSize = 8;
    // - - - - - - - - - - - - - - - - - - -

    // We use 8 bits per fragment
    ExplicitBitVect* fp = new ExplicitBitVect(mNChunks * chunkSize);
    // ...
    RDKit::MorganFingerprints::BitInfoMap atomsSettingBits;
    // We need only fragments, not the fingerprint.
    RDKit::SparseIntVect<boost::uint32_t>* bitVector = RDKit::MorganFingerprints::getFingerprint(
            *mol, mRadius, 0, 0, false, false, true, false, &atomsSettingBits);
    delete bitVector;
    // Find all fragments.
    for (RDKit::MorganFingerprints::BitInfoMap::iterator iter = atomsSettingBits.begin();
            iter != atomsSettingBits.end(); ++iter)
    {
        boost::uint32_t position = iter->first % mNChunks;
        VectorInfo& infoVector = iter->second;
        // Iterate over fragments info.
        for (VectorInfo::iterator itemIter = infoVector.begin(); itemIter != infoVector.end(); ++itemIter)
        {
            if (itemIter->second == 0)
            {
                // Ignore zero size.
                continue;
            }
            // Assemble fragments.
            RDKit::PATH_TYPE bondIndexes = RDKit::findAtomEnvironmentOfRadiusN(*mol,
                    itemIter->second, itemIter->first);
            std::vector<int> atomsToUse;
            for (RDKit::PATH_TYPE::iterator bondIter = bondIndexes.begin(); bondIter != bondIndexes.end();
                    ++bondIter)
            {
                // Store indexes.
                atomsToUse.push_back(mol->getBondWithIdx(*bondIter)->getBeginAtomIdx());
                atomsToUse.push_back(mol->getBondWithIdx(*bondIter)->getEndAtomIdx());
            }
            std::sort(atomsToUse.begin(), atomsToUse.end());
            atomsToUse.erase(std::unique(atomsToUse.begin(), atomsToUse.end()), atomsToUse.end());
            // Get fragment as smile.
            std::string fragmentSmile = RDKit::MolFragmentToSmiles(*mol, atomsToUse, &bondIndexes);

            if (cache.count(fragmentSmile) == 0)
            {
                // Convert fragment into molecule.
                RDKit::ROMol *fragment = RDKit::SmilesToMol(fragmentSmile);
                // Compute some statistic and set respective bits to one.
                cache[fragmentSmile] = boost::shared_ptr<ExplicitBitVect>(new ExplicitBitVect(chunkSize));
                // - - - - - VECTORFP SETTINGS - - - - -


                // Best parameterization from Article:                
                unsigned int nHBDon_Lipinski = RDKit::Descriptors::calcNumHBD(*fragment);
                unsigned int nN = std::count(fragmentSmile.begin(), fragmentSmile.end(), 'N');

                storeValueBinaryEncoding(nHBDon_Lipinski, 4, *cache[fragmentSmile].get(), 0);
                storeValueBinaryEncoding(nN, 3, *cache[fragmentSmile].get(), 4);
                
                // Set final bit for fragment presence, so we do not get empty fingerprint.
                cache[fragmentSmile].get()->setBit(chunkSize - 1);

                // - - - - - - - - - - - - - - - - - - -
                delete fragment;
                fragment = 0;
            }
            // Get value from cache and store into fingerprint.
            int basePosition = position * chunkSize;
            // Copy and store on right position.
            std::vector<int> bitsOn;
            cache[fragmentSmile].get()->getOnBits(bitsOn);
            for (std::vector<int>::iterator iter = bitsOn.begin(); iter != bitsOn.end(); ++iter) {
                fp->setBit(basePosition + *iter);
            }
        }
    }
    return fp;    
}

