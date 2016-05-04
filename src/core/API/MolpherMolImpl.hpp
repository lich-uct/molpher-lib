/* 
 * File:   MolpherMolImpl.hpp
 * Author: sichom
 *
 * Created on January 27, 2016, 3:47 PM
 */

#ifndef MOLPHERMOLIMPL_HPP
#define	MOLPHERMOLIMPL_HPP

#include <memory>
#include <vector>

#include "core/data_structs/MolpherMolData.hpp"
#include "data_structs/MolpherMol.hpp"
#include "data_structs/ExplorationTree.hpp"

class MolpherMol::MolpherMolImpl {
    
    friend class MolpherMol;
    
private:
    std::shared_ptr<ExplorationTree> tree;
    MolpherMolData data;

public:
    MolpherMolImpl(const std::string& SMILES);
    MolpherMolImpl(const MolpherMolData& data);
    MolpherMolImpl(const MolpherMolImpl& other);
    MolpherMolImpl();
    
    std::unique_ptr<MolpherMolImpl> copy() const;

    void setSMILES(const std::string& smiles);
};

#endif	/* MOLPHERMOLIMPL_HPP */

