/**
    \file morphing_functions.hpp

    Include `morphing_functions.hpp` to use miscellaneous morphing functions
    that wrap some high level functionality of Molpher.

    As of yet the file contains the definition of a single function
    which uses a simplified implementation of the Morphing algorithm
    used in Molpher v1.01.
*/

#ifndef TESTING_HPP
#define	TESTING_HPP

#include <string>

/**
 * Runs the original morphing algorithm
 * using the BasicPathFinder implementation.
 *
 * @param storagePath name of the folder where the results are saved
 * @param jobFile path to a job file or XML template
 * @param threadCnt max. number of threads to create
 */
void run_path_finder(
    const std::string &storagePath
    , const std::string &jobFile
    , int threadCnt
);

#endif	/* TESTING_HPP */

