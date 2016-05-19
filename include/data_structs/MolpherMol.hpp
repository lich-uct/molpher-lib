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

#ifndef MOLPHERMOL_HPP
#define	MOLPHERMOL_HPP

#include <string>
#include <memory>
#include <set>

#include "selectors/chemoper_selectors.h"

class ExplorationTree; // forward declaration to resolve circular dependency

class MolpherMol {
    
private:
    class MolpherMolImpl;
    std::unique_ptr<MolpherMolImpl> pimpl;

public:  
    MolpherMol();
    MolpherMol(const std::string& smiles, const std::string& formula, const std::string& parentSmile,
                const unsigned& oper, const double& dist, const double& distToClosestDecoy,
                const double& weight, const double& sascore);
    MolpherMol(const std::string& smiles);
    MolpherMol(const MolpherMol& other);
    ~MolpherMol();
    
    MolpherMol& operator=(const MolpherMol&);
    std::shared_ptr<MolpherMol> copy() const;
    
    // getters
    const std::string& getSMILES() const;
    double getDistToTarget() const;
    std::shared_ptr<ExplorationTree> getTree();
    const std::string& getParentSMILES() const;
    const std::set<std::string>& getDescendants() const;
    const std::set<std::string>& getHistoricDescendants() const;
    unsigned int getItersWithoutDistImprovement() const;
    double getSAScore() const;
    double getMolecularWeight() const;
    const std::string& getFormula() const;
    int getParentOper() const;
    
    // setters
    void setOwner(std::shared_ptr<ExplorationTree> tree);
    void setSMILES(const std::string&);
    void setParentSMILES(const std::string&);
    void setDistToTarget(double dist);
    void setSAScore(double score);
    void setItersWithoutDistImprovement(unsigned int count);
    void increaseItersWithoutDistImprovement();
    void decreaseItersWithoutDistImprovement();
    void addToDescendants(const std::string& smiles);
    void removeFromDescendants(const std::string& smiles);
    void setDescendants(const std::set<std::string>&);
    void addToHistoricDescendants(const std::string& smiles);
    void removeFromHistoricDescendants(const std::string& smiles);
    void setHistoricDescendants(const std::set<std::string>&);
    
    bool isValid() const;
    bool isBoundToTree() const;
    void removeFromTree();
};

#endif	/* MOLPHERMOL_HPP */

