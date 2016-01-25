/*
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author: Petyr
 * Created on 18.10.2013
 */

#include <string.h>
#include <iostream>
#include <fstream>

#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include "boost/archive/archive_exception.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/xml_parser.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/BadFileException.h>

#include "inout.h"
#include "iteration_serializer.hpp"
#include "selectors/fingerprint_selectors.h"
#include "selectors/simcoeff_selectors.h"

namespace molpher {
namespace iteration {

enum FileType {
    XML_FILE,
    XML_TEMPLATE_FILE,
    SNP_FILE,
    UNKNOWN_FILE
};

/**
 * Save snapshot in snp file format.
 * @param file
 * @param snp
 * @return
 */
bool saveSnp(const std::string &file, const IterationSnapshot &snp) {
    std::ofstream outStream;
    outStream.open(file.c_str(), std::ios_base::out | std::ios_base::trunc);
    if (!outStream.good()) {
        return false;
    }

    bool failed = false;
    try {
        boost::archive::text_oarchive oArchive(outStream);
        oArchive << snp;
    } catch (boost::archive::archive_exception &exc) {
        failed = true;
        SynchCout(std::string(exc.what()));
    }
    outStream.close();
    return !failed;
}

/**
 * Save snapshot as xml.
 * @param file
 * @param snp
 * @return
 */
bool saveXml(const std::string &file, const IterationSnapshot &snp) {
    std::ofstream outStream;
    outStream.open(file.c_str(), std::ios_base::out | std::ios_base::trunc);
    if (!outStream.good()) {
        return false;
    }

    bool failed = false;
    try {
        boost::archive::xml_oarchive oArchive(outStream);
        oArchive << BOOST_SERIALIZATION_NVP(snp);
    } catch (boost::archive::archive_exception &exc) {
        failed = true;
        SynchCout(std::string(exc.what()));
    }
    outStream.close();
    return !failed;
}

/**
 * Load snapshot from snp file.
 * @param file
 * @param snp
 * @return
 */
bool loadSnp(const std::string &file, IterationSnapshot &snp) {
    std::ifstream inStream;
    inStream.open(file.c_str());
    if (!inStream.good()) {
        return false;
    }

    try {
        boost::archive::text_iarchive iArchive(inStream);
        iArchive >> snp;
    } catch (boost::archive::archive_exception &exc) {
        SynchCout(std::string(exc.what()));
        inStream.close();
        return false;
    }

    inStream.close();
    return true;
}

/**
 * Load snapshot from xml file.
 * @param file
 * @param snp
 * @return
 */
bool loadXml(const std::string &file, IterationSnapshot &snp) {
    std::ifstream inStream;
    inStream.open(file.c_str());
    if (!inStream.good()) {
        return false;
    }

    try {
        boost::archive::xml_iarchive iArchive(inStream);
        iArchive >> BOOST_SERIALIZATION_NVP(snp);
    } catch (boost::archive::archive_exception &exc) {
        SynchCout(std::string(exc.what()));
        inStream.close();
        return false;
    }

    inStream.close();
    return true;
}

int toInt(const std::string& str) {
    std::stringstream ss;
    ss << str;
    int result;
    ss >> result;
    return result;
}

/**
 * Create molecule from given smile. The smile may change as it's
 * RDKit smile
 * @param inSmile
 * @return
 */
MolpherMolecule createMoleculeFromSmile(const std::string& inSmile) {
    RDKit::RWMol* mol = RDKit::SmilesToMol(inSmile);
    try {
        RDKit::MolOps::Kekulize(*mol);
    } catch (const ValueErrorException &exc) {
        SynchCout("Cannot kekulize input molecule.");
    }

    std::string smile(RDKit::MolToSmiles(*mol));
    std::string formula(RDKit::Descriptors::calcMolFormula(*mol));

    SynchCout("Parse molecule " + inSmile + " >> " + smile);

    return MolpherMolecule(smile, formula);
}

/**
 * Based on given template fetch data into given snapshot. Values
 * that are not specified in the template are ignored.
 * @param is
 * @param snp
 */
void loadXmlTemplate(std::istream &is, IterationSnapshot &snp) {
    boost::property_tree::ptree pt;
    // read xml
    boost::property_tree::read_xml(is, pt);

    BOOST_FOREACH( boost::property_tree::ptree::value_type const& v, pt.get_child("iteration") ) {

        if (v.first == "source") {
            snp.source = createMoleculeFromSmile(v.second.get<std::string>("smile"));
        } else if (v.first == "target") {
            snp.target = createMoleculeFromSmile(v.second.get<std::string>("smile"));
        } else if (v.first == "fingerprint") {
            snp.fingerprintSelector = FingerprintParse(v.second.data());
        } else if (v.first == "similarity") {
            snp.simCoeffSelector = SimCoeffParse(v.second.data());
        } else if (v.first == "param") {
            BOOST_FOREACH( boost::property_tree::ptree::value_type const& v,
                    v.second ) {
                // parameters ..
                if (v.first == "syntetizedFeasibility") {
                    snp.params.useSyntetizedFeasibility =
                            v.second.data() == "1" || v.second.data() == "true";
                } else if (v.first == "acceptMin") {
                    snp.params.cntCandidatesToKeep = toInt(v.second.data());
                } else if (v.first == "acceptMax") {
                    snp.params.cntCandidatesToKeepMax = toInt(v.second.data());
                } else if (v.first == "farProduce") {
                    snp.params.cntMorphs = toInt(v.second.data());
                } else if (v.first == "closeProduce") {
                    snp.params.cntMorphsInDepth = toInt(v.second.data());
                } else if (v.first == "farCloseThreashold") {
                    std::stringstream ss;
                    ss << v.second.data();
                    ss >> snp.params.distToTargetDepthSwitch;
                } else if (v.first == "maxMorhpsTotal") {
                    snp.params.cntMaxMorphs = toInt(v.second.data());
                } else if (v.first == "nonProducingSurvive") {
                    snp.params.itThreshold = toInt(v.second.data());
                } else if (v.first == "iterMax") {
                    snp.params.cntIterations = toInt(v.second.data());
                } else if (v.first == "maxTimeMinutes") {
                    snp.params.timeMaxSeconds = toInt(v.second.data()) * 60;
                } else if (v.first == "weightMin") {
                    snp.params.minAcceptableMolecularWeight = toInt(v.second.data());
                } else if (v.first == "weightMax") {
                    snp.params.maxAcceptableMolecularWeight = toInt(v.second.data());
                }
            }
        } else {
            // unexpected token
        }
    }

    // printout some statistics .. first prepare the report
    std::stringstream ss;
    ss << "The new iteration has been created from template: " << std::endl;
    ss << "\tsource: " << snp.source.smile << std::endl;
    ss << "\ttarget: " << snp.target.smile << std::endl;
    // and print ..
    SynchCout(ss.str());
}

bool loadXmlTemplate(const std::string &file, IterationSnapshot &snp) {
    std::ifstream inStream;
    inStream.open(file.c_str(), std::ios_base::in);
    if (!inStream.good()) {
        return false;
    }
    loadXmlTemplate(inStream, snp);
    inStream.close();
    return true;
}

/**
 * Return type of file based on extension.
 * @param file
 * @return File type.
 */
FileType fileType(const std::string &file) {
    // we start with the more specific ..
    if (boost::algorithm::ends_with(file, "-template.xml")) {
        return XML_TEMPLATE_FILE;
    } else if (boost::algorithm::ends_with(file, ".snp")) {
        return SNP_FILE;
    } else if (boost::algorithm::ends_with(file, ".xml")) {
        return XML_FILE;
    } else {
        return UNKNOWN_FILE;
    }
}

bool IterationSerializer::save(const std::string &file, const IterationSnapshot &snp) const {
    switch(fileType(file)) {
        default:
        case UNKNOWN_FILE:
        case SNP_FILE:
            return saveSnp(file, snp);
        case XML_FILE:
            return saveXml(file, snp);
        case XML_TEMPLATE_FILE:
            // we do not support save in xml template format
            return false;
    }
}

bool IterationSerializer::load(const std::string &file, IterationSnapshot &snp) const {
    switch(fileType(file)) {
        case SNP_FILE:
            return loadSnp(file, snp);
        case XML_FILE:
            return loadXml(file, snp);
        case XML_TEMPLATE_FILE:
            return loadXmlTemplate(file, snp);
        case UNKNOWN_FILE:
        default:
            return false;
    }
}

} }