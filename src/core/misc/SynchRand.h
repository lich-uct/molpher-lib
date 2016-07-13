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

#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>

#include <tbb/spin_mutex.h>
#include <tbb/enumerable_thread_specific.h>

// [0..3] - 0 most precise, 3 fastest
#ifndef SYNCHRAND_SPEED_OVER_PRECISION
#define SYNCHRAND_SPEED_OVER_PRECISION 3
#endif
 
/**
 * Synchronized random number generator.
 * 
 * Class should be used as a singleton.
 */
class SynchRand
{
public:
    /**
     * Generates uniformly distributed random number from a given range.
     * @param min
     * @param max
     * @return 
     */    
    static int GetRandomNumber(int min, int max);
    /**
     * Generates uniformly distributed non-negative random number.
     * @param max
     * @return 
     */
    static int GetRandomNumber(int max);

    static void SetSeed(unsigned seed);
private:
    SynchRand();
    SynchRand(const SynchRand &other);
    SynchRand &operator=(const SynchRand &other);    
private:
#if SYNCHRAND_SPEED_OVER_PRECISION == 0
    typedef boost::random::mt19937 Engine;
    typedef boost::random::uniform_int_distribution<int> Distribution;
#elif SYNCHRAND_SPEED_OVER_PRECISION == 1
    typedef boost::random::mt11213b Engine;
    typedef boost::random::uniform_int_distribution<int> Distribution;
#elif SYNCHRAND_SPEED_OVER_PRECISION == 2
    typedef boost::random::mt19937 Engine;
    typedef boost::random::uniform_smallint<int> Distribution;
#elif SYNCHRAND_SPEED_OVER_PRECISION == 3
    typedef boost::random::mt11213b Engine;
    typedef boost::random::uniform_smallint<int> Distribution;
#else
    #error SynchRand: Precision not specified.
#endif    
    typedef tbb::spin_mutex Guard;    
    typedef tbb::enumerable_thread_specific<Engine> Engines;    
private:    
    /**
     * Singleton instance.
     */    
    static SynchRand instance;
    /**
     * Mutex.
     */
    Guard mSeedGuard;
    
    boost::random::mt19937 mSeedEngine;
    
    boost::random::uniform_int_distribution<int> mSeedDistribution;
    
    Engines mEngines;
};


