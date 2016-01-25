
#include "callbacks/EraseSubtreeCallback.hpp"

EraseSubtreeCallback::EraseSubtreeCallback(PathFinderContext& context) : TraverseCallback(context) {
    // no action
}


void EraseSubtreeCallback::processMorph(MolpherMol& morph) {
    std::deque<std::string> toErase;
    toErase.push_back(morph.getSMILES());

    while (!toErase.empty()) {
        std::string current = toErase.front();
        toErase.pop_front();

        PathFinderContext::CandidateMap::accessor ac;
        context->candidates.find(ac, current);
        assert(!ac.empty());

        std::set<std::string>::const_iterator it;
        for (it = ac->second.descendants.begin();
                it != ac->second.descendants.end(); it++) {
            toErase.push_back(*it);
        }

        context->prunedDuringThisIter.push_back(current);
        context->candidates.erase(ac);
    }
}
