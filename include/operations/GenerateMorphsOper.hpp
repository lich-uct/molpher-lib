/* 
 * File:   PutativeExtend.hpp
 * Author: sichom
 *
 * Created on October 19, 2015, 8:39 AM
 */

#ifndef PUTATIVEEXTEND_HPP
#define	PUTATIVEEXTEND_HPP

#include <memory>

#include "TreeOperation.hpp"

class GenerateMorphsOper : public TreeOperation {
    
public:
    class GenerateMorphsOperImpl;
    GenerateMorphsOper(std::shared_ptr<ExplorationTree> expTree);
    GenerateMorphsOper();
    virtual void operator()();
    
private:
    std::shared_ptr<GenerateMorphsOperImpl> pimpl;

};

#endif	/* PUTATIVEEXTEND_HPP */

