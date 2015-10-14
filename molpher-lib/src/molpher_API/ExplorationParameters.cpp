
#include "molpher_API/ExplorationParameters.hpp"

ExplorationParameters::ExplorationParameters() {
    // no action
}

IterationSnapshot ExplorationParameters::createIterationSnapshot() {
    return iterSnapshot;
}
