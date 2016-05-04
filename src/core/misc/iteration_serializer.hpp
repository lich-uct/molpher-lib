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

#pragma once

#include "data_structs/ExplorationData.hpp"
#include "core/data_structs/ExplorationDataImpl.hpp"

namespace molpher {
namespace iteration {
    
//MolpherMolecule createMoleculeFromSmile(const std::string& inSmile);

/**
 * Provide functionality that enable save/load iteration snapshot 
 * into file.
 */
class IterationSerializer
{    
public:
    /**
     * Save snapshot into file.
     * @param file Name of file into which save.
     * @param snp Snapshot to save.
     */
    static bool save(const std::string &file, const ExplorationData::ExplorationDataImpl &data);
    /**
     * Load snapshot from file. Based on given file extension decide
     * which file to load.
     * @param file Name of file from that load data.
     * @param snp Snapshot into which load data.
     * @return
     */
    static bool load(const std::string &file, ExplorationData::ExplorationDataImpl &data);
};

} }