/*
 Copyright (c) 2012 Petr Koupy

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

#include <ctime>

#include "SynchRand.h"

SynchRand::SynchRand()
{
#if SYNCHRAND_DETERMINISTIC == 1
    mSeedEngine.seed(0);
#else
    // TODO - on MPI cluster, node ID should be included into seed
    mSeedEngine.seed(static_cast<unsigned int>(std::time(NULL)));
#endif
}

SynchRand SynchRand::instance;

int SynchRand::GetRandomNumber(int min, int max)
{
    Distribution distribution(min, max);

    // Retrieve engine from Thread Local Storage.
    bool alreadyInitializedEngine = false;
    Engines::reference engine = instance.mEngines.local(alreadyInitializedEngine);

    if (!alreadyInitializedEngine) {
        /* This slow path will run at most once for each created worker thread.
         According to TBB documentation, worker thread count should be in most
         of the cases same as the number of physical computational units in
         the system. Because worker threads are reused for TBB tasks, it could
         be even said that slow path will be executed at most once for each
         physical computational unit. */
        Guard::scoped_lock lock(instance.mSeedGuard);
        engine.seed(instance.mSeedDistribution(instance.mSeedEngine));
    }

    return distribution(engine);
}

int SynchRand::GetRandomNumber(int max)
{
    return SynchRand::GetRandomNumber(0, max);
}
