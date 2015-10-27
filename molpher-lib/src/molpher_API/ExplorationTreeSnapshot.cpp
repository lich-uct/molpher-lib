
#include <stdexcept>
#include <iostream>

#include "molpher_API/ExplorationTreeSnapshot.hpp"
#include "inout.h"

ExplorationTreeSnapshot::ExplorationTreeSnapshot(IterationSnapshot& iterSnapshot) {
    if (iterSnapshot.IsValid()) {
        this->iterSnapshot = iterSnapshot;
    } else {
        throw std::runtime_error("Supplied IterationSnapshot is invalid."); // TODO wrap this in SWIG
    }
}

ExplorationTreeSnapshot* ExplorationTreeSnapshot::load(const std::string& filename) {
    IterationSnapshot iterSnap;
    if (ReadSnapshotFromFile(filename, iterSnap) && iterSnap.IsValid()) {
        std::cout << "Snapshot successfully created from: " << filename << std::endl;
    } else {
        throw std::runtime_error(std::string("Supplied IterationSnapshot is invalid. Tree snapshot cannot be opened:").append(filename)); // TODO wrap this in SWIG
    }
    return new ExplorationTreeSnapshot(iterSnap);
}

void ExplorationTreeSnapshot::save(const std::string& filename) {
    WriteSnapshotToFile(filename, iterSnapshot);
}
