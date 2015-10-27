
#include <tbb/task_scheduler_init.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_sort.h>

#include "molpher_API/operations/SortMorphsOper.hpp"

SortMoprhsOper::SortMoprhsOper(ExplorationTree& expTree) : TreeOperation(expTree) {
    // no action
}

bool SortMoprhsOper::CompareMorphs::operator()(
        const MolpherMolecule &a, const MolpherMolecule &b) const {
    /* Morphs are rated according to their proximity to the connecting line
     between their closest decoy and the target (i.e. sum of both distances is
     minimal on the connecting line between decoy and target). When sums for
     both morphs are equal, it is possible (but not necessary) that both
     morphs lie on the same connecting line. In that case, morphs are
     rated only according to their proximity to the target. Such comparison
     should allow convergence to the target even in the late stages of the
     algorithm when majority of morphs lie on the connecting line between
     decoy closest to the target and the target itself. */

    double aSum = a.distToTarget + a.distToClosestDecoy;
    double bSum = b.distToTarget + b.distToClosestDecoy;

    bool approximatelyEqual = (
            fabs(aSum - bSum) <= (32 * DBL_EPSILON * fmax(fabs(aSum), fabs(bSum))));

    if (approximatelyEqual) {
        return a.distToTarget < b.distToTarget;
    } else {
        return aSum < bSum;
    }

    /**
     * Just sort based on distance to decoy .. or target. We prefer those
     * with target.
     */

    /* Experimental evaluation
    if (a.nextDecoy == -1 && b.nextDecoy == -1) {
        // both go for target, so take just their distance
        return a.distToTarget < b.distToTarget;
    } else if (a.nextDecoy == -1) {
        // a goes for target while b not -> a is greater
        return true;
    } else if (b.nextDecoy == -1) {
        // b goes for target while a not -> b is greater
        return false;
    } else {
        // decide based on nextDecoy
        if (a.nextDecoy == b.nextDecoy) {
            // go for same decoy use distance to closes decoy
            return a.distToClosestDecoy < b.distToClosestDecoy;
        } else {
            // greater is better
            return a.nextDecoy > b.nextDecoy;
        }
    }*/
}

void SortMoprhsOper::operator()() {
    tbb::task_scheduler_init scheduler;
    if (threadCnt > 0) {
        scheduler.terminate();
        scheduler.initialize(threadCnt);
    }
    
    ExplorationTree::MoleculeVector& morphs = fetchGeneratedMorphs();
    CompareMorphs compareMorphs;
    tbb::parallel_sort(morphs.begin(), morphs.end(), compareMorphs);
}

